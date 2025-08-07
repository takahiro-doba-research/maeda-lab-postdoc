from .coordinate import Coordinate


class Output:

    def __init__(
            self, name=None, coordinate=None, E=None, spin=None,
            gradient=None, hessian=None, normal_mode=None, freq=None,
            E_zpv=None, E_E_zpv=None,
            E_tr=None, E_rot=None, E_vib=None, H_corr=None, H=None,
            S_el=None, S_tr=None, S_rot=None, S_vib=None, G_corr=None, G=None,
            E_zpv_r=None, E_E_zpv_r=None,
            E_tr_r=None, E_rot_r=None, E_vib_r=None, H_corr_r=None, H_r=None,
            S_el_r=None, S_tr_r=None, S_rot_r=None, S_vib_r=None,
            G_corr_r=None, G_r=None,
            connection=None):

        self.name = name
        self.coordinate = coordinate

        if not isinstance(self.coordinate, Coordinate):
            self.coordinate = Coordinate()
        
        self.E = E                    # a
        self.spin = spin
        self.gradient = gradient
        self.hessian = hessian
        self.normal_mode = normal_mode
        self.freq = freq

        # Thermochemistry
        self.E_zpv = E_zpv            # b
        self.E_E_zpv = E_E_zpv        # c = a + b
        self.E_tr = E_tr              # d
        self.E_rot = E_rot            # e
        self.E_vib = E_vib            # f
        self.H_corr = H_corr          # g = d + e + f + boltzmann_const * temp
        self.H = H                    # h = a + g
        self.S_el = S_el              # i
        self.S_tr = S_tr              # j
        self.S_rot = S_rot            # k
        self.S_vib = S_vib            # l
        self.G_corr = G_corr          # m = g - (i + j + k + l) * temp
        self.G = G                    # n = a + m

        # Thermochemistry after the replacement of small eigenvalues.
        self.E_zpv_r = E_zpv_r
        self.E_E_zpv_r = E_E_zpv_r
        self.E_tr_r = E_tr_r
        self.E_rot_r = E_rot_r
        self.E_vib_r = E_vib_r
        self.H_corr_r = H_corr_r
        self.H_r = H_r
        self.S_el_r = S_el_r
        self.S_tr_r = S_tr_r
        self.S_rot_r = S_rot_r
        self.S_vib_r = S_vib_r
        self.G_corr_r = G_corr_r
        self.G_r = G_r

        self.connection = connection

    def show(self):
        self.coordinate.show(self.name)

    def save(self, path, header=None, footer=None, ignore_notes=False):
        
        if header == None:
            header = ["# STO-3G\n", "\n", "title\n", "\n", "0 1\n"]
        if footer == None:
            footer = ["\n"]

        self.coordinate.save(path, header, footer, ignore_notes)
