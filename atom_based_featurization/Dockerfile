FROM python:3.11.14
ENV PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update && \
    apt-get install -y less locales vim && \
    localedef -f UTF-8 -i ja_JP ja_JP.UTF-8
ENV LANG=ja_JP.UTF-8 \
    LANGUAGE=ja_JP:ja \
    LC_ALL=ja_JP.UTF-8 \
    TZ=Asia/Tokyo \
    TERM=xterm

RUN wget https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz && \
    tar -xf xtb-6.7.1-linux-x86_64.tar.xz && \
    mkdir -p /opt/xtb && \
    mv xtb-dist/* /opt/xtb/ && \
    rm -rf xtb-dist xtb-6.7.1-linux-x86_64.tar.xz

ENV PATH="/opt/xtb/bin:$PATH"

COPY requirements.txt .
RUN pip install --upgrade pip setuptools && \
    pip install --no-cache-dir -r requirements.txt

COPY grrmlib ./grrmlib
RUN pip install -e ./grrmlib