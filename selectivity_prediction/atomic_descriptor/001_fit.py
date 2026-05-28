import numpy as np
import polars as pl
from sklearn.linear_model import Ridge
from sklearn.metrics import root_mean_squared_error
from sklearn.preprocessing import StandardScaler

file = "001_fit"

df = pl.read_parquet("dataset.parquet")

column_group = "pyridone"
column_fix = "backbone"
columns_feat = [col for col in df.columns if col.startswith(("pyridone_", "backbone_"))]
column_target = "beta_av_logit"

df_train_val_preds = []
df_test_preds = []
df_coefs = []

"""
Outer cross-validation
"""
groups_test = df[column_group].unique(maintain_order=True)

for group_test in groups_test:
    print(f"group_test: {group_test}")
    
    df_train_val = df.filter(pl.col(column_group) != group_test)
    df_test = df.filter(pl.col(column_group) == group_test)
    
    """
    Recursive feature elimination (RFE)
    """
    columns_rfe = columns_feat.copy()
    
    while len(columns_rfe) > 4:
        n_features = len(columns_rfe)
        print(f"n_features: {n_features}")
        
        X_train_val = df_train_val.select(columns_rfe).to_numpy()
        y_train_val = df_train_val.select(column_target).to_numpy()
        X_test = df_test.select(columns_rfe).to_numpy()
        y_test = df_test.select(column_target).to_numpy()
        
        ss_X = StandardScaler()
        ss_y = StandardScaler()
        
        X_train_val_ = ss_X.fit_transform(X_train_val)
        y_train_val_ = ss_y.fit_transform(y_train_val)
        X_test_ = ss_X.transform(X_test)
        y_test_ = ss_y.transform(y_test)
        
        """
        Hyperparameter optimization
        """
        alpha = 100
        
        for i in range(1, -3, -1):
            rmses = {}
            
            for j in range(-9, 10):
                y_val_preds_ = []
                y_val_trues_ = []
                
                """
                Inner cross-validation
                """
                groups_val = df_train_val[column_group].unique(maintain_order=True)
                
                for group_val in groups_val:
                    df_train = df_train_val.filter(pl.col(column_group) != group_val)
                    df_val = df_train_val.filter(pl.col(column_group) == group_val)
                    
                    X_train = df_train.select(columns_rfe).to_numpy()
                    y_train = df_train.select(column_target).to_numpy()
                    X_val = df_val.select(columns_rfe).to_numpy()
                    y_val = df_val.select(column_target).to_numpy()
                    
                    X_train_ = ss_X.transform(X_train)
                    y_train_ = ss_y.transform(y_train)
                    X_val_ = ss_X.transform(X_val)
                    y_val_ = ss_y.transform(y_val)
                    
                    model = Ridge(alpha + j * 10**i)
                    model.fit(X_train_, y_train_)
                    
                    y_val_pred_ = model.predict(X_val_)
                    y_val_preds_.append(y_val_pred_)
                    y_val_trues_.append(y_val_)
                
                y_val_preds_ = np.concatenate(y_val_preds_)
                y_val_trues_ = np.concatenate(y_val_trues_)
                rmses[j] = root_mean_squared_error(y_val_trues_, y_val_preds_)
            
            j_min = min(rmses, key=rmses.get)
            alpha += j_min * 10**i
        
        """
        Refit with the best hyperparameter
        """
        model = Ridge(alpha)
        model.fit(X_train_val_, y_train_val_)
        
        """
        Store data
        """
        y_train_val_pred_ = model.predict(X_train_val_)
        y_train_val_pred = ss_y.inverse_transform(y_train_val_pred_)
        df_train_val_pred = pl.DataFrame(
            {
                f"{column_group}_test": group_test,
                "n_features": n_features,
                column_group: df_train_val[column_group],
                column_fix: df_train_val[column_fix],
                column_target: df_train_val[column_target],
                f"{column_target}_pred": y_train_val_pred.ravel()
            }
        )
        df_train_val_preds.append(df_train_val_pred)
        
        y_test_pred_ = model.predict(X_test_)
        y_test_pred = ss_y.inverse_transform(y_test_pred_)
        df_test_pred = pl.DataFrame(
            {
                f"{column_group}_test": group_test,
                "n_features": n_features,
                column_fix: df_test[column_fix],
                column_target: df_test[column_target],
                f"{column_target}_pred": y_test_pred.ravel()
            }
        )
        df_test_preds.append(df_test_pred)
        
        coef = model.coef_
        df_coef = (
            pl.DataFrame(coef, schema=columns_rfe)
            .with_columns(
                pl.lit(group_test).alias(f"{column_group}_test"),
                pl.lit(n_features).alias("n_features"),
                pl.lit(alpha).alias("alpha"),
            )
            .select([f"{column_group}_test", "n_features", "alpha"] + columns_rfe)
        )
        df_coefs.append(df_coef)
        
        """
        Feature elimination
        """
        column_eliminate = columns_rfe[np.argmin(np.abs(coef))]
        columns_rfe = [column for column in columns_rfe if column != column_eliminate]

df_train_val_preds = pl.concat(df_train_val_preds, how="vertical_relaxed")
df_test_preds = pl.concat(df_test_preds, how="vertical_relaxed")
df_coefs = pl.concat(df_coefs, how="diagonal_relaxed")

df_train_val_preds.write_parquet(f"{file}_train_val_preds.parquet")
df_test_preds.write_parquet(f"{file}_test_preds.parquet")
df_coefs.write_parquet(f"{file}_coefs.parquet")