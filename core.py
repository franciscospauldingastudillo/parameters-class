class parameters:
    def __init__(self,**params):
        # Spectral resolution used in RFM experiments
        self.nu0 = params.get('nu0', 10)  # cm-1
        self.nu1 = params.get('nu1', 1500)  # cm-1
        self.dnu = params.get('dnu', 0.1)  # cm-1
        self.exp = params.get('exp', 'earth')
        self.band = params.get('band', 'wv-broadband')
        self.runtype = params.get('runtype', 'cooling')
        self.cpdef = params.get('cp', 29012 / 28.964)  # J/kg/K (default RFM value of specific heat)
        self.nsday = params.get('nsday', 86400)  # seconds per Earth-day
        self.TEMREL = params.get('TEMREL', 0)  # ground-air temp diff.

        # spectral range calculations
        self.nus = np.arange(self.nu0, self.nu1 + self.dnu, self.dnu)
        self.nnus = len(self.nus)
        self.update_spectral_range()

        # Thermodynamic parameters
        self.ps = params.get('ps', 1e5)
        self.Gamma = params.get('Gamma', 6.5e-3)
        self.z = params.get('z', np.arange(0, 5e4, 1))

        # Saturation Vapor Pressure (Water Vapor)
        self.Ttrip = params.get('Ttrip', 273.16)  # K 
        self.ptrip = params.get('ptrip', 611.65)  # Pa
        self.E0v   = params.get('E0v', 2.3740e6)  # J/kg 
        self.ggr   = params.get('ggr', 9.81)      # m/s^2, gravity
        self.rgasa = params.get('rgasa', 287.04)  # J/kg/K, specific gas constant of dry air
        self.rgasv = params.get('rgasv', 461.5)   # J/kg/K, specific gas constant of water vapor
        self.cva   = params.get('cva', 719)       # J/kg/K, isovolumetric specific heat of dry air
        self.cvv   = params.get('cvv', 1418)      # J/kg/K, isovolumetric specific heat of water vapor
        self.cvl   = params.get('cvl', 4119)      # J/kg/K, isovolumetric specific heat of liquid water
        self.cpa   = params.get('cpa', self.cva + self.rgasa)   # isobaric specific heat of dry air
        self.cpv   = params.get('cpv', self.cvv + self.rgasv)   # isobaric specific heat of water vapor
        self.eps   = params.get('eps', self.rgasa / self.rgasv) # ratio of specific gas constants
        self.L     = params.get('L', 2.5e6)       # enthalpy of vaporization of water
        self.E0s   = params.get('E0s', np.nan)    # no ice phase
        self.cvs   = params.get('cvs', np.nan)    # no ice phase
           
        # Earth-like composition parameters
        self.Runi  = params.get('Runi',8.314)               # universal gas constant (J/mol/K)
        self.MN2   = params.get('MN2',2*14.01*1e-3)         # molecular weight of N2 (kg/mol)
        self.MO2   = params.get('MO2',2*15.999*1e-3)        # molecular weight of O2 (kg/mol)
        self.MH2O  = params.get('MH2O',18.015*1e-3)         # molecular weight of H2O (kg/mol)
        self.MCO2  = params.get('MCO2',44.009*1e-3)         # molecular weight of CO2 (kg/mol)
        self.RN2   = params.get('RN2',self.Runi/self.MN2)   # specific gas constant of N2 (J/kg/K)
        self.RO2   = params.get('RO2',self.Runi/self.MO2)   # specific gas constant of O2 (J/kg/K)
        self.RH2O  = params.get('RH2O',self.Runi/self.MH2O) # specific gas constant of H2O (J/kg/K)
        self.RCO2  = params.get('RCO2',self.Runi/self.MCO2) # specific gas constant of CO2 (J/kg/K)
        # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Type=JANAFG&Plot=on
        self.cpmolN2  = params.get('cpmolN2',29.10)         # molar specific heat of N2 (J/mol/K)
        self.cpmolO2  = params.get('cpmolO2',29.30)         # molar specific heat of O2 (J/mol/K)
        self.cpmolCO2 = params.get('cpmolCO2',37.22)        # molar specific heat of CO2 (J/mol/K)
        self.cpmolH2O = params.get('cpmolH2O',33.58)        # molar specific heat of H2O (J/mol/K)
        
        self.cpN2  = params.get('cpN2',self.cpmolN2/self.MN2)      # specific heat of N2 (J/kg/K)
        self.cpO2  = params.get('cpO2',self.cpmolO2/self.MO2)      # specific heat of O2 (J/kg/K)
        self.cpCO2 = params.get('cpCO2',self.cpmolCO2/self.MCO2)   # specific heat of CO2 (J/kg/K)
        self.cpH2O = params.get('cpH2O',self.cpmolH2O/self.MH2O)   # specific heat of H2O (J/kg/K)
        
        # Earth-like molar mixing ratios
        self.xCO2  = params.get('xCO2',395e-6)         # molar mixing ratio of CO2 (mol/mol) ~ 395 ppmv
        self.xN2   = params.get('xN2',0.78)            # molar mixing ratio of N2 (mol/mol)  ~ 78% (dry-air approx)
        self.xO2   = params.get('xO2',0.22)            # molar mixing ratio of O2 (mol/mol)  ~ 22% (dry-air approx)
        
        # redundancies
        self.Rd   = 287
        self.Rv   = self.RH2O

        # SSM1D
        self.use_ssm1d_ctm = False
        self.update_ssm1d_params()

        # Initialize case string
        self.case = ""

    def generate_case(self,**params):
        exp          = params.get('exp','unknown')
        self.Ts      = params.get('Ts', 300)
        self.Tmid    = params.get('Tmid', 250)
        self.Ttrp    = params.get('Ttrp', 200)
        self.RHs     = params.get('RHs', 0.75)
        self.RHmid   = params.get('RHmid', 0.54)
        self.RHtrp   = params.get('RHtrp', 0.75)
        self.uniform = params.get('uniform',1)
        self.update_alpha()
        # Extract gases dynamically (gas1, gas2, ..., gasN)
        gases = [params[key] for key in sorted(params) if key.startswith('gas')]
        # Define valid CIA pairs as tuples for safety with multi-character molecules
        self.valid_ciapairs = params.get('valid_ciapairs', [('N2', 'N2'), ('N2', 'CH4'), ('N2', 'H2'), ('CH4', 'CH4')])
        # Filter CIA pairs to include only those with gases in the `gases` list
        ciapairs = [
            f"{mol1}{mol2}" for (mol1, mol2) in self.valid_ciapairs
            if mol1 in gases and mol2 in gases
        ]
        # Format the CIA part of the case name
        cia_str = "-CIA-" + "-".join(ciapairs) if ciapairs else ""
        # Construct the case name
        self.case = '-'.join([
            exp,
            *gases,
            f"Ts{int(self.Ts)}",
            f"Tmid{int(self.Tmid)}",
            f"Ttrp{int(self.Ttrp)}",
            f"RHs{int(self.RHs * 100)}",
            f"RHmid{int(self.RHmid * 100)}",
            f"RHtrp{int(self.RHtrp * 100)}",
            str(self.uniform)
        ]) + cia_str
        print('generate_case: ',self.case)
        
    def update_spectral_range(self):
        """ Update spectral range and band indices. """        
        if self.band == 'wv-rot':
            self.i0 = np.squeeze(np.where(self.nus == self.nu0))
            self.i1 = np.squeeze(np.where(np.abs(self.nus - 1000) == np.min(np.abs(self.nus - 1000))))
        elif self.band == 'wv-rot-right':
            self.i0 = np.squeeze(np.where(np.abs(self.nus - 150) == np.min(np.abs(self.nus - 150))))
            self.i1 = np.squeeze(np.where(np.abs(self.nus - 1000) == np.min(np.abs(self.nus - 1000))))
        elif self.band == 'wv-vib-rot':
            self.i0 = np.squeeze(np.where(np.abs(self.nus - 1000) == np.min(np.abs(self.nus - 1000))))
            self.i1 = np.squeeze(np.where(np.abs(self.nus - self.nu1) == np.min(np.abs(self.nus - self.nu1))))
        elif self.band == 'wv-broadband':
            self.i0 = np.squeeze(np.where(np.abs(self.nus - self.nu0) == np.min(np.abs(self.nus - self.nu0))))
            self.i1 = np.squeeze(np.where(np.abs(self.nus - self.nu1) == np.min(np.abs(self.nus - self.nu1))))
        elif self.band == 'co2':
            self.i0 = np.squeeze(np.where(self.nus == self.nu0))
            self.i1 = np.squeeze(np.where(np.abs(self.nus - self.nu1) == np.min(np.abs(self.nus - self.nu1))))
        elif self.band == 'ch4':
            self.i0 = np.squeeze(np.where(self.nus == self.nu0))
            self.i1 = np.squeeze(np.where(np.abs(self.nus - self.nu1) == np.min(np.abs(self.nus - self.nu1))))
        else:
            raise ValueError(f"Invalid band type: {self.band}. Expected 'earth'.")

    def update_band(self, new_band):
        """ Update the band and recalculate spectral range values. """
        self.band = new_band
        self.update_spectral_range()

    def update_alpha(self):
        # Equation (13) from Spaulding-Astudillo and Mitchell (2025)
        self.alpha_lt = -np.log(self.RHtrp / self.RHmid)/(self.Ttrp - self.Tmid) ** 2
        self.alpha_gt = -np.log(self.RHs / self.RHmid)/(self.Ts - self.Tmid) ** 2

    def update_use_ssm1d_ctm(self, new_use_ssm1d_ctm):
        """ Update SSM1D continuum transmission model. """
        self.use_ssm1d_ctm = new_use_ssm1d_ctm
        self.update_ssm1d_params()

    def update_ssm1d_params(self):
        """ Update SSM1D regression parameters. """
        if not self.use_ssm1d_ctm: # without continuum
            print('switching to SSM1D-NOCTM')
            self.D    = 1.5 # two-steam diffusivity coefficient
            self.pref = 50000
            self.Tref = 260 
            self.pinf  = 2.5e11 # JF20 assume 2.5e11
            # rotation band
            self.nurot = 150 # cm-1
            self.lrot0 = 56.06 # cm-1
            self.krot  = 139.57 # m2/kg
            # line-only component
            self.Prot_lo = 375 # cm-1
            self.Arot_lo = 0.15
            # continuum-only component
            self.Prot_co = 1000-150 # cm-1
            self.Arot_co = 0
            # vibration band
            self.nuvib = 1500 # cm-1
            self.lvib0 = 39.53 # cm-1
            self.kvib  = 11.78 # m2/kg
            # line-only component
            Pvib_lo = 0
            Avib_lo = 0
            # continuum-only component
            Pvib_co = 0
            Avib_co = 0
        else: # with continuum
            print('switching to SSM1D-CTM')
            self.D    = 1.5 # two-steam diffusivity coefficient
            self.pref = 50000
            self.Tref = 260 
            self.pinf  = 1.0e11 # JF20 assume 2.5e11
            # rotation band
            self.nurot = 150 # cm-1
            self.lrot0 = 74.36 # cm-1
            self.krot  = 78.77 # m2/kg
            # line-only component
            self.Prot_lo = 375 # cm-1
            self.Arot_lo = 0.15
            # continuum-only component
            self.Prot_co = 1000-150 # cm-1
            self.Arot_co = 0.55
            # vibration band
            self.nuvib = 1500 # cm-1
            self.lvib0 = 39.53 # cm-1
            self.kvib  = 11.78 # m2/kg
            # line-only component
            Pvib_lo = 0
            Avib_lo = 0
            # continuum-only component
            Pvib_co = 0
            Avib_co = 0
