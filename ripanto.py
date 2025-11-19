from matplotlib import lines
import numpy as np
import os
import re

class RipartoSpinelli:
    def __init__(
        self,JK,valenz,el,bc,nomi
    ):
        self.JK = int(JK)
        self.VALENZ = np.asarray(valenz, dtype=float)
        self.EL = np.asarray(el, dtype=float)
        self.BC = np.asarray(bc, dtype=float)
        self.NOMI = list(nomi)

        n = self.JK
        if not (len(self.VALENZ) == len(self.EL) == len(self.BC) == len(self.NOMI) == n):
            raise ValueError(
                f"All arrays must have length JK={n} "
                f"(got VALENZ={len(self.VALENZ)}, EL={len(self.EL)}, "
                f"BC={len(self.BC)}, NOMI={len(self.NOMI)})"
            )

    @classmethod
    def with_defaults(cls):
        JK = 22
        VALENZ = [3,2,3,2,2,3,4,3,2,3,2,4,0,2.5,2,0,0,0,0,0,0,0]
        EL = [13.0,26.0,26.0,12.0,25.0,25.0,14.0,24.0,30.0,23.0,
              28.0,22.0,0.0,26.0,27.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        BC = [3.449,9.54,9.54,5.375,-3.73,-3.73,4.149,3.635,5.689,
              -0.3824,10.3,-3.438,0.0,9.54,2.50,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        NOMI = [
            '  Al3+ ', '  Fe2+ ', '  Fe3+ ', '  Mg2+ ', '  Mn2+ ',
            '  Mn3+ ', '  Si4+ ', '  Cr3+ ', '  Zn2+ ', '  V3+  ', '  Ni2+ ',
            '  Ti4+ ', '  Vac  ', '  Fe25 ', '  Co2+ ', '  In3+ ', '  Ga3+ ',
            '  Sn4+ ', '  Mo4+ ', '  Ag2+ ', '  Cu2+ ', '  xxx  '
        ]
        return cls(JK, VALENZ, EL, BC, NOMI)

    @staticmethod
    def extract_title(input_file, number = 89):
        """
        Extracts a formatted title string from the input file name.
        Example: 'SP-3G12.IN' with number=89 -> 'Spinello 3g12 raff. pri89'
        """
        fname = os.path.basename(input_file)
        fname_no_ext = os.path.splitext(fname)[0]
        match = re.search(r'(\d+[A-Z]+\d*)', fname_no_ext)
        code = match.group(0) if match else ""
        return f"Spinello {code.lower()} raff. pri{number}"
    
    @staticmethod
    def prepare_input_file(input_file):
        """
        Prepares the input file to ensure it has the correct format
        to be read by the FCN function.
        """
        with open(input_file, "r", encoding="utf-8") as f:
            lines = f.readlines()

        new_lines = lines[:1] + lines[17:]

        with open(input_file, "w", encoding="utf-8") as f:
            f.writelines(new_lines)


    def FCN(self, X, IFLAG, input_file, dis_file):
        if IFLAG == 1:
            # --- Step 1: Read input file ---
            with open(input_file, 'r') as fin:
                lines = fin.readlines()
            # Find IXN line and parse parameters
            # ixn_line = [l for l in lines if 'IXN' in l]
            # ixn_idx = lines.index(ixn_line[0])
            self.IXN, self.TOL, self.INOR, self.TTT, self.PPP = map(float, lines[1].split())
            self.DELTAT = self.TTT - 25
            self.DELTAP = self.PPP - 0.0001
            
            self.XXX = np.copy(X)
            self.title = self.extract_title(input_file)

    # names = tuple(n.strip() for n in self.NOMI if n.strip())
    # cation_names = []
    # for l in lines:
    #     s = l.lstrip()                  
    #     if s.startswith(names):         # starts with ANY of the names
    #         cation_names = s.split()    
    #         break
            self.cation_names = lines[3].split()
            self.K9 = len(self.cation_names)

            # --- Step 2: Read DIS file ---
            with open(dis_file, 'r') as fdis:
                dis_lines = fdis.readlines()
            self.dis_names = dis_lines[0].split()
            
           #dis_name_to_index = {name: i for i, name in enumerate(dis_names)}
            self.K7 = len(self.dis_names)

            self.DDT = np.array([float(x) for x in dis_lines[1].split()])
            self.CCDT = np.array([float(x) for x in dis_lines[2].split()])
            self.dDTdp = np.array([float(x) for x in dis_lines[3].split()])
            self.dDTdt = np.array([float(x) for x in dis_lines[4].split()])
            self.DDM = np.array([float(x) for x in dis_lines[5].split()])
            self.CCDM = np.array([float(x) for x in dis_lines[6].split()])
            self.dDMdp = np.array([float(x) for x in dis_lines[7].split()])
            self.dDMdt = np.array([float(x) for x in dis_lines[8].split()])
            self.COST1T, self.COST2T, self.HP1, self.HP2, self.HP3 = map(float, dis_lines[9].split()[:5])
            self.SIGMAA, self.SIGMAU, self.TOLL, self.INORR, self.SIGCH = map(float, dis_lines[10].split()[:5])
            self.SST, self.SSM, self.SSQ = map(float, dis_lines[11].split()[:3])
            self.IMINU = int(dis_lines[12].split()[0])
            print(' iminu=', self.IMINU)

                    # --- Step 3: Main calculation logic (translated) ---
        # Fill missing dDTdt and dDMdt using Hazen & Prewitt
            for j in range(self.K7):
                if self.dDTdt[j] < 1e-8 and self.VALENZ[j] > 0.1:
                    # 4.0 is the coordination number
                    self.ALFAT = self.COST1T * (self.COST2T - self.VALENZ[j]/4.0) * 1e-6
                    if self.ALFAT <= 0.0:
                        self.ALFAT = 0.0
                    self.dDTdt[j] = self.ALFAT * self.DDT[j] # apply correction to the distance

                if self.dDMdt[j] < 1e-8 and self.VALENZ[j] > 0.1:
                    self.ALFAM = self.COST1T * (self.COST2T - self.VALENZ[j]/6.0) * 1e-6
                    self.dDMdt[j] = self.ALFAM * self.DDM[j] # apply correction to the distance

                if self.dDTdp[j] < 1e-8 and self.VALENZ[j] > 0.1:
                    self.BETAT = 37.0 * (self.DDT[j]**3 / self.VALENZ[j]) * 1e-5
                    self.dDTdp[j] = self.BETAT * self.DDT[j] # correction for pressure

                if self.dDMdp[j] < 1e-8 and self.VALENZ[j] > 0.1:
                    self.BETAM = 37.0 * (self.DDM[j]**3 / self.VALENZ[j]) * 1e-5
                    self.dDMdp[j] = self.BETAM * self.DDM[j] # correction for pressure

            self.DDT = self.DDT - self.DELTAP * self.dDTdp + self.DELTAT * self.dDTdt # recomputed distances
            self.DDM = self.DDM - self.DELTAP * self.dDMdp + self.DELTAT * self.dDMdt
            self.IOP1 = 0
            self.IOP2 = 0
            self.VALE = np.zeros(self.K9)
            self.ELET = np.zeros(self.K9)
            self.DT = np.zeros(self.K9)
            self.DM = np.zeros(self.K9)
            self.CDT = np.zeros(self.K9)
            self.CDM = np.zeros(self.K9)
            for i in range(self.K9):
                ih = 0
                for j in range(self.K7):
                    if self.cation_names[i] == self.dis_names[j]:
                        ih += 1
                        self.VALE[i] = self.VALENZ[j]
                        self.ELET[i] = self.EL[j] if self.IXN == 1 else self.BC[j]
                        self.DT[i] = self.DDT[j]
                        self.DM[i] = self.DDM[j]
                        self.CDT[i] = self.CCDT[j]
                        self.CDM[i] = self.CCDM[j]
                        if self.cation_names[i] == 'Fe2+':
                            self.IOP1 = i
                        if self.cation_names[i] == 'Fe3+':
                            self.IOP2 = i
                if ih == 0:
                    raise RuntimeError(f" CATIONE {self.cation_names[i]} NON TROVATO IN FILE DISTANZE")

            # Read AT and ERR arrays
           # at_idx = [i for i, l in enumerate(lines) if l.strip().startswith('Al3+')][0] + 1
            self.AT = np.array([float(x) for x in lines[4].split()]) # 'bulk' occupancies
            print("AT occupancies:", *self.AT[:self.K9])
            self.ERR = np.array([float(x) for x in lines[5].split()])
            print("ERR occupancies:", *self.ERR[:self.K9])
            self.KAT = np.sum(self.AT[:self.K9] > 0.02)
            print(f" Number of significant cations KAT={self.KAT}")

            self.ET, self.EM, self.A, self.U, self.SIGET, self.SIGEM, self.SIGA, self.SIGU = map(float, lines[7].split())

            if self.SIGMAA > 1e-7:
                self.SIGA = self.SIGMAA

            if self.SIGMAU > 1e-7:
                self.SIGU = self.SIGMAU

            if self.SIGCH > 1e-5:
                self.ERR[:self.K9] = self.SIGCH * self.AT[:self.K9]

            if self.TOLL > 1e-4:
                self.TOL = self.TOLL

            if self.INORR == 1:
                self.INOR = 1

            self.ST, self.SM, self.SQ = map(float, lines[9].split())
        # Modify cation distances for T and P
        #ZDT = DDT.copy()
        #ZDM = DDM.copy()
            if self.SST > 1e-7:
                self.ST = self.SST
            if self.SSM > 1e-7:
                self.SM = self.SSM
            if self.SSQ > 1e-7:
                self.SQ = self.SSQ
            if self.SIGET < 1e-8:
                self.KAT -= 1
            if self.SIGEM < 1e-8:
                self.KAT -= 1

            SUMCAT = np.sum(self.AT[:13])
            if abs(SUMCAT - 3.0) >= 0.001:
                print(f"I CATIONI CHIUDONO A {SUMCAT}")

            CARICA = np.sum(self.AT[:self.K9] * self.VALE[:self.K9])
            if abs(CARICA - 8.0) >= 0.001:
                print(f"LE CARICHE CHIUDONO A {CARICA}")

            # using neutron data (not X ray/ electronic data),
            # so the electron counts are not adjusted
            if self.IXN == 2:
                self.ET0 = self.ET
                self.EM0 = self.EM
            else:
                self.ELCHIM = np.sum(self.AT[:self.K9] * self.ELET[:self.K9])
                self.ELSTR = self.ET + 2.0 * self.EM
                FATT = (self.ELCHIM - self.ELSTR) / self.ELCHIM
                print(f"   e- chim., e- strutt., fattore(CH.-ST.) {self.ELCHIM:.3f} {self.ELSTR:.3f} {FATT:.3f}")
                self.ET0 = self.ET
                self.EM0 = self.EM
                if self.INOR != 0:
                    self.ET += FATT * self.ET
                    self.EM += FATT * self.EM
            self.DTO = self.A * np.sqrt(3.0 * (0.125 - self.U) ** 2)
            self.DMO = self.A * np.sqrt((0.5 - self.U) ** 2 + 2.0 * (0.25 - self.U) ** 2)
            self.DER1 = np.sqrt(3.0 * (0.125 - self.U) ** 2) # partial derivative of DTO w.r.t. A
            self.DER2 = (self.A * np.sqrt(3) * (-0.25 + 2 * self.U)) / (2.0 * np.sqrt(3.0 * (0.125 - self.U) ** 2)) # partial derivative of DMO w.r.t. U
            self.SIGTO = np.sqrt(self.DER1 ** 2 * self.SIGA ** 2 + self.DER2 ** 2 * self.SIGU ** 2)
            #Combined uncertainty on DTO, propagated from uncertainties in A and U using standard error propagation.       
            self.DER3 = np.sqrt((0.5 - self.U) ** 2 + 2.0 * (0.25 - self.U) ** 2) # partial derivative of DMO w.r.t. A
            self.DER4 = (self.A * (6.0 * self.U - 2.0)) / (2.0 * np.sqrt((0.5 - self.U) ** 2 + 2.0 * (0.25 - self.U) ** 2)) # partial derivative of DMO w.r.t. U
            self.SIGMO = np.sqrt(self.DER3 ** 2 * self.SIGA ** 2 + self.DER4 ** 2 * self.SIGU ** 2)
            # Combined uncertainty on DMO, propagated from uncertainties in A and U.
        # DISTANCE FILE WINS OVER INPUT FILE
        # If the DIS file provides a non-negligible value for SIGMAA or SIGMAU, 
        # use it to override the input file's error on a and u.
        
            
            F = self.calc_res(X, IFLAG, input_file)
            return F
            
        else:
            F = self.calc_res(X, IFLAG, input_file)
            return F


    def calc_res(self, X, IFLAG, input_file):
            # --- Section: Calculate residuals (lines 285-420) ---
            SC = np.zeros(8) # this will store the main residuals
        # differences between calculated and expected values
        # for various constants

        # eM (electrons in M site)
        # X is the array of cation site occupancies being optimized.
            EEM = np.sum(X[self.K9:self.K9*2] * self.ELET[:self.K9])
            SC[0] = EEM - 2.0 * self.EM

        # Charge of M (QQM) and occupancy of M (OCM)
            QQM = np.sum(X[self.K9:self.K9*2] * self.VALE[:self.K9])
            OCM = np.sum(X[self.K9:self.K9*2])
            SC[1] = OCM - 2.0

        # eT (electrons in T site)
            EET = np.sum(X[:self.K9] * self.ELET[:self.K9])
            SC[2] = EET - self.ET

        # Charge of T (QQT) and occupancy of T (OCT)
            QQT = np.sum(X[:self.K9] * self.VALE[:self.K9])
            OCT = np.sum(X[:self.K9])
            SC[3] = OCT - 1.0

        # Calculate DTCAL (T-O distance, including M-site contribution)
            DTCAL = np.sum(X[:self.K9] * self.DT[:self.K9]) # average T-O distance 
            DTCAL += np.sum(X[self.K9:self.K9*2] * self.CDT[:self.K9]) # average T-O distance with M-site corrective contribution 

        # Calculate DMCAL (M-O distance, including Fe2.5+ logic)
            DMCAL = 0.0
            FE25M = 0.0

        # Fe2.5+ Logic (Marshall & Dollase electron hopping)
        #Purpose: If both Fe2+ and 
        # Fe3+ are present in the M site 
        # (octahedral site) and their occupancies 
        # are similar (within a tolerance TOL), 
        # they are combined into Fe2.5+ (representing
        #  electron hopping). How:
        #   If their difference is small, 
        #combine all into Fe2.5+ and set Fe2+ and Fe3+ occupancies to zero.
        # If not, combine as much as possible into Fe2.5+, 
        # leaving the remainder as either Fe2+ or Fe3+.
            if self.TOL <= 90.0 and self.IOP1 > 0 and self.IOP2 > 0 and X[self.IOP1+self.K9] > 0.1 and X[self.IOP2+self.K9] > 0.1:
                FEM = X[self.IOP1+self.K9] + X[self.IOP2+self.K9]
                SALV1 = X[self.IOP1+self.K9]
                SALV2 = X[self.IOP2+self.K9]
                QMAX = max(X[self.IOP1+self.K9], X[self.IOP2+self.K9])
                if abs(X[self.IOP1+self.K9] - X[self.IOP2+self.K9]) <= (self.TOL * QMAX):
                    FE25M = X[self.IOP1+self.K9] + X[self.IOP2+self.K9]
                    X[self.IOP1+self.K9] = 0.0
                    X[self.IOP2+self.K9] = 0.0
                elif X[self.IOP1+self.K9] < X[self.IOP2+self.K9]:
                    FE25M = 2.0 * (X[self.IOP1+self.K9] + self.TOL * X[self.IOP1+self.K9])
                    X[self.IOP1+self.K9] = 0.0
                    X[self.IOP2+self.K9] = FEM - FE25M
                elif X[self.IOP1+self.K9] > X[self.IOP2+self.K9]:
                    FE25M = 2.0 * (X[self.IOP2+self.K9] + self.TOL * X[self.IOP2+self.K9])
                    X[self.IOP2+self.K9] = 0.0
                    X[self.IOP1+self.K9] = FEM - FE25M

            DMCAL = np.sum(X[self.K9:self.K9*2] * self.DM[:self.K9])

            if FE25M >= 0.001:
                X[self.IOP1+self.K9] = SALV1 #they were set to zero above
                X[self.IOP2+self.K9] = SALV2

            DMCAL += FE25M * self.DDM[14] #TO DO: save the index and use it instead of hardcoding 14

            DMCAL += np.sum(X[:self.K9] * self.CDM[:self.K9]) # average M-O distance with T-site corrective contribution
            DMCAL /= 2.0 # average M-O distance (since M there are two M sites per formula unit in spinel)

        # Calculate AO and UC
            RADI = 33.0 * DMCAL ** 2 - 8.0 * DTCAL ** 2
            SC[4] = 0.0

            if RADI >= 0.0:
                AO = (np.sqrt(RADI) + 5.0 * DTCAL) * 0.4198911
                # The factor 0.4198911 converts this sum to the cell edge length (Ångström).
                SC[4] = AO - self.A

            if abs(DTCAL - DMCAL) < 1e-6:
                UC = 0.2625 # set it to a default value
            # which is typical for normal spinels
            # R = (DMCAL/DTCAL)**2 = 1
            else:
                R = (DMCAL ** 2) / (DTCAL ** 2)
                RQ = (33.0 / 16.0) * R - 0.5
                if RQ < 0.0:
                    print("RADICE PER IL CALCOLO DI U NEGATIVA!")
                    RQ = 0.0
                UC = (0.75 * R - 2.0 + np.sqrt(RQ)) / (6.0 * R - 6.0)
            SC[5] = UC - self.U
        # UC : model-predicted oxygen positional parameter

        # Calculate DTCAL and DMCAL residuals
        # --- Translation of RIPANTO.FOR lines 421-698 ---
        # Calculate ATT, CARIC, SC[7], SC[8], and F
        # Calculate ATT (sum of T and M site for each cation)
            ATT = X[:self.K9] + X[self.K9:self.K9*2]
        # aka total occupancy of each cation

        # Calculate total charge
            CARIC = np.sum(self.VALE[:self.K9] * ATT)
            SC[6] = CARIC - 8.0

        # Calculate chi-square residual for cation distribution
            SC[7] = 0.0
            for K in range(self.K9):
                if self.ERR[K] < 1e-9:
                    self.ERR[K] = 0.0010
                SC[7] += ((self.AT[K] - ATT[K]) / self.ERR[K]) ** 2 # AT is given as input from XRAY data
            # while ATT is calculated from the occupancies X that we want to optimize

# Calculate objective function F
            if self.SIGEM < 1e-8:
                SCAR1 = 0.0
            else:
                SCAR1 = (0.5 * SC[0] / self.SIGEM) ** 2

            if self.SIGET < 1e-8:
                SCAR3 = 0.0
            else:
                SCAR3 = (SC[2] / self.SIGET) ** 2

            F = (
            SCAR1
            + (0.5 * SC[1] / self.SM) ** 2
            + SCAR3
            + (SC[3] / self.ST) ** 2
            + (SC[4] / self.SIGA) ** 2
            + (SC[5] / self.SIGU) ** 2
            + (SC[6] / self.SQ) ** 2
            + SC[7]
            )
            F = F / (7.0 + self.KAT)
# You can implement file writing and formatted output here if needed.
          # --- Step 5: Faithful Fortran-style file output (IFLAG == 3 case) ---
            if IFLAG == 3:
            # 2. ARCHI.AGG append at end (simulate BACKSPACE)
                SUMAT = np.sum(self.AT)
                SUMATT = np.sum(ATT)
                foutname = os.path.basename(input_file)
                foutname_no_ext = os.path.splitext(foutname)[0]
                
                with open(foutname_no_ext + ".OUT", "a") as f17:
                    f17.write("*******************************************************\n")
                    f17.write(f"{self.title}\n")
                    f17.write(f"F(Xi)={F:10.3e} Norm. e-: {self.INOR}  TOLL% HOP={self.TOL:6.2f}  T°C {self.TTT:5.0f}  PS(GPa) {self.PPP:10.4f}\n")
                    if self.IXN == 1:
                        f17.write(f"         e-T    e-M    A0        u        T-O     M-O       e-X     e-chim\n")
                        f17.write(f"OBS {self.ET0:7.3f} {self.EM0:7.3f} {self.A:9.5f} {self.U:10.6f} {self.DTO:8.4f} {self.DMO:8.4f} {self.ET+2*self.EM:9.3f} {self.ELCHIM:9.3f}\n")
                    elif self.IXN == 2:
                        f17.write(f"  BCT    BCM    A0        u        T-O     M-O\n")
                        f17.write(f"OBS {self.ET0:7.3f} {self.EM0:7.3f} {self.A:9.5f} {self.U:10.6f} {self.DTO:8.4f} {self.DMO:8.4f}\n")
                    f17.write(f"SIG {self.SIGET:7.3f} {self.SIGEM:7.3f} {self.SIGA:9.5f} {self.SIGU:10.6f} {self.SIGTO:8.4f} {self.SIGMO:8.4f} car(T)   car(M)\n")
                    f17.write(f"CAL {EET:7.3f} {EEM/2.0:7.3f} {AO:9.5f} {UC:10.6f} {DTCAL:8.4f} {DMCAL:8.4f} {QQT:8.4f} {QQM/2.0:8.4f}\n")
                    f17.write(f"SC  {EET-self.ET:7.3f} {EEM/2.0-self.EM:7.3f} {AO-self.A:9.5f} {UC-self.U:10.6f} {DTCAL-self.DTO:8.4f} {DMCAL-self.DMO:8.4f}   CARICA   {SC[6]:8.4f}\n")
                    ERT = np.sqrt((X[:self.K9] / ATT[:self.K9]) * self.ERR[:self.K9] ** 2)
                    ERM = np.sqrt((X[self.K9:self.K9*2] / ATT[:self.K9]) * self.ERR[:self.K9] ** 2)
                    # T-site occupancies
                    f17.write("    " + " ".join(f"{name:7s}" for name in self.cation_names[:self.K9]) + " SOMMA\n")
                    f17.write("T " + " ".join(f"{val:7.4f}" for val in X[:self.K9]) + f" {OCT:7.4f}\n")
                    # T-site errors
                    f17.write("  " + " ".join(f"{val:7.4f}" for val in ERT) + "\n")
                    # M-site occupancies
                    f17.write("M " + " ".join(f"{val:7.4f}" for val in X[self.K9:self.K9*2]) + f" {OCM:7.4f}\n")
                    # M-site errors
                    f17.write("  " + " ".join(f"{val:7.4f}" for val in ERM) + "\n")
                    if FE25M >= 0.001:
                        f17.write(f" Fe2.5+: {FE25M:7.4f}  TOL (Marshall e Dollase) = {self.TOL:7.4f}\n")
                    if self.TOL > 98:
                        f17.write(" Non si vuole calcolare Fe2.5+ (TOL = 99)\n")
                        # Observed occupancies
                    f17.write("      " + " ".join(f"{name:7s}" for name in self.cation_names[:self.K9]) + " SOMMA\n")
                    f17.write(" obs " + " ".join(f"{val:7.4f}" for val in self.AT[:self.K9]) + f" {SUMAT:7.4f}\n")
                    # Calculated occupancies
                    f17.write(" cal " + " ".join(f"{val:7.4f}" for val in ATT[:self.K9]) + f" {SUMATT:7.4f}\n")
                    # Errors
                    f17.write(" sig " + " ".join(f"{val:7.4f}" for val in self.ERR[:self.K9]) + "\n")
                    kk = 7 + self.KAT
                    if self.IXN == 1:
                        # Electron calculation summary
                        f17.write(
                        f" Chi(i) quadro;{kk:2d} residui, compresi {self.KAT:2d} cationi >0.02\n"
                        "  eT       OccT     eM       OccM     A0       U        neut     chim\n"
                        f"{SCAR3:9.2e} {(SC[3]/self.ST)**2:9.2e} {SCAR1:9.2e} {(0.5*SC[1]/self.SM)**2:9.2e} "
                        f"{(SC[4]/self.SIGA)**2:9.2e} {(SC[5]/self.SIGU)**2:9.2e} {(SC[6]/self.SQ)**2:9.2e} {SC[7]:9.2e}\n"
                    )
                    elif self.IXN == 2:
                # Bond compressibility summary
                        f17.write(
                            f" Chi(i) quadro;{kk:2d} residui, compresi {self.KAT:2d} cationi >0.02\n"
                            " BCT       OccT    BCM       OccM     A0       U        neut     chim\n"
                            f"{SCAR3:9.2e} {(SC[3]/self.ST)**2:9.2e} {SCAR1:9.2e} {(0.5*SC[1]/self.SM)**2:9.2e} "
                            f"{(SC[4]/self.SIGA)**2:9.2e} {(SC[5]/self.SIGU)**2:9.2e} {(SC[6]/self.SQ)**2:9.2e} {SC[7]:9.2e}\n"
                            )

                    f17.write("*******************************************************\n")


                with open("ARCHI.AGG", "a+") as f18:
                    f18.write(f"{self.title}\n")
                    f18.write(" ".join([f"{val:7.4f}" for val in self.AT[:self.K9]]) + "\n")
                    f18.write(" ".join([f"{val:7.4f}" for val in self.ERR[:self.K9]]) + "\n")
                    f18.write("T " + " ".join([f"{val:7.4f}" for val in X[:self.K9]]) + "\n")
                    f18.write("M " + " ".join([f"{val:7.4f}" for val in X[self.K9:self.K9*2]]) + "\n")
                    f18.write("  e-T    e-M    A0     u    SIGMA eT     eM     a0     u\n")
                    f18.write(f"{self.ET0:7.3f} {self.EM0:7.3f} {self.A:8.5f} {self.U:8.5f} {self.SIGET:7.3f} {self.SIGEM:7.3f} {self.SIGA:8.5f} {self.SIGU:8.5f}\n")
                    f18.write("    " + " ".join(f"{name:7s}" for name in self.cation_names[:self.K9]) + " SOMMA\n")
                    f18.write(" " + " ".join(f"{val:7.4f}" for val in self.AT[:self.K9]) + "\n")
                    f18.write(" " + " ".join(f"{val:7.4f}" for val in self.ERR[:self.K9]) + "\n")
                    f18.write("T" + "".join(f"{val:7.4f}" for val in X[:self.K9]) + "\n")
                    f18.write("M" + "".join(f"{val:7.4f}" for val in X[self.K9:self.K9*2]) + "\n")
                    f18.write("  e-T    e-M    A0     u    SIGMA eT     eM     a0     u\n")
                    f18.write(f" {self.ET0:7.3f} {self.EM0:7.3f} {self.A:8.5f} {self.U:8.5f} {self.SIGET:7.3f} {self.SIGEM:7.3f} {self.SIGA:8.5f} {self.SIGU:8.5f}\n")


                # ...existing code...
                FR = [0.0] * (2 * self.K7)   # Output array for calculated occupancies
                FR0 = [0.0] * (2 * self.K7)  # Output array for initial occupancies

                for i in range(self.K7):
                    found = False
                    for j in range(self.K9):
                        # Compare names after stripping spaces
                        if self.dis_names[i].strip() == self.cation_names[j].strip():
                            # T-site (first K7)
                            FR[i] = X[j]
                            FR0[i] = self.XXX[j]
                            # M-site (second K7)
                            FR[i + self.K7] = X[j + self.K9]
                            FR0[i + self.K7] = self.XXX[j + self.K9]
                            found = True
                            break
                  #  if not found:
                       # print(f"{self.dis_names[i]} NON TROVATA CORRISPONDENZA !")
            


                with open("FORXLS.OBS", "r+") as f16:
                    f16lines = f16.readlines()
                    KRIGHE = len(f16lines)
                    if KRIGHE == 0:
                        f16.write("A0     u      T-O     M-O    e-T    e-M" +
                            "".join([f"T{n:7s}" for n in self.NOMI[:15]]) +
                            "".join([f"M{n:7s}" for n in self.NOMI[:15]]) + "\n")
                        f16.write(f" {self.A:8.5f} {self.U:7.4f} {self.DTO:7.4f} {self.DMO:7.4f} {self.ET0:7.3f} {self.EM0:7.3f}" +
                            "".join([f"{val:8.4f}" for val in FR0[:30]]) + "\n")

                # 4. FORXLS.CAL overwrite
                with open("FORXLS.CAL", "r+") as f23:
                    f23lines = f23.readlines()
                    KRIGHE = len(f23lines)
                    if KRIGHE == 0:
                        f23.write("FCN     A0     u      T-O     M-O    e-T    e-M" +
                            "".join([f"T{n:7s}" for n in self.NOMI[:15]]) +
                            "".join([f"M{n:7s}" for n in self.NOMI[:15]]) + "\n")
                        f23.write(f"{F:7.2f} {AO:8.5f} {UC:7.4f} {DTCAL:7.4f} {DMCAL:7.4f} {EET:7.3f} {EEM/2.0:7.3f}" +
                            "".join([f"{val:8.4f}" for val in FR[:30]]) + "\n")

                if self.IMINU == 1:
                    self.aggiorna(X, input_file)
            return F

        # --- Step 4: Output results (minimal, for now just return F) ---
        # The full output logic (writing files) can be implemented as needed.
        # For now, just return a dummy F for demonstration.
    def aggiorna(self, X, input_file='MINUIT.IN'):
        """
        Aggiorna il file MINUIT.IN (o altro specificato) con i dati correnti dell'oggetto.
        """
        with open(input_file, 'r+') as f:
            lines = f.readlines()
            TIT = lines[0].strip() if lines else ""
            f.seek(0)
            f.write(f"{TIT}\n")
            TOT = np.sum(self.AT[:12])

            AMAX = np.zeros(self.K9)
            for i in range(self.K9):
                if self.AT[i] >= 1.0:
                    AMAX[i] = self.AT[i] + 0.004 * self.AT[i]
                elif 1.0 > self.AT[i] > 0.4:
                    AMAX[i] = self.AT[i] + 0.005 * self.AT[i]
                elif 0.4 >= self.AT[i] >= 0.1:
                    AMAX[i] = self.AT[i] + 0.008 * self.AT[i]
                elif 0.1 > self.AT[i] >= 0.05:
                    AMAX[i] = self.AT[i] + 0.010 * self.AT[i]
                elif 0.05 > self.AT[i] >= 0.02:
                    AMAX[i] = self.AT[i] + 0.025 * self.AT[i]
                elif 0.02 > self.AT[i] >= 0.005:
                    AMAX[i] = self.AT[i] + 0.050 * self.AT[i]
                elif self.AT[i] < 0.005:
                    AMAX[i] = self.AT[i] + 0.250 * self.AT[i]

            # Scrivi cationi T-site
            for i in range(self.K9):
                f.write(f"{i+1:3.0f}       {self.cation_names[i]:10s}{X[i]:10.4f}{X[i]/60.:10.4f}{0.0:10.4f}{AMAX[i]:10.4f}\n")
            # Scrivi cationi M-site
            for i in range(self.K9):
                f.write(f"{i+1+self.K9:3.0f}       {self.cation_names[i]:10s}{X[i+self.K9]:10.4f}{X[i+self.K9]/60.:10.4f}{0.0:10.4f}{AMAX[i]:10.4f}\n")

            f.write("\n")
            f.write("IXN,0/(1,2,3) PER HOP., NORM. e-,TEMP, PRES \n")
            f.write(f"{int(self.IXN):2d}{self.TOL:7.3f}{int(self.INOR):3d}{self.TTT:8.1f}{self.PPP:8.4f}\n")
            f.write("Cationi in CAT(I), cariche in CAR(I), errori ASSOLUTI in ERR(I)\n")
            # Scrivi nomi cationi
            for i in range(0, self.K9, 10):
                f.write("".join(f"{self.cation_names[j]:7s}" for j in range(i, min(i+10, self.K9))) + "\n")
            # Scrivi AT
            for i in range(0, self.K9, 13):
                f.write(" ".join(f"{self.AT[j]:7.4f}" for j in range(i, min(i+13, self.K9))) + "\n")
            # Scrivi ERR
            for i in range(0, self.K9, 13):
                f.write(" ".join(f"{self.ERR[j]:7.4f}" for j in range(i, min(i+13, self.K9))) + "\n")
            f.write("  eT     eM     a0      u      Sigma:eT   eM    a0    u\n")
            f.write(f"{self.ET:7.3f}{self.EM:7.3f}{self.A:9.5f}{self.U:10.6f}{self.SIGET:7.3f}{self.SIGEM:7.3f}{self.SIGA:9.5f}{self.SIGU:9.5f}\n")
            f.write(f"SigmaT  SigmaM  SigmaCar in ST, SM, SQ\n{self.ST:9.5f}{self.SM:9.5f}{self.SQ:9.5f}\n")
            f.write("\n")
            f.write("PRINTOUT        0.\nSIMPLEX     20000.\nMIGRAD      45000.   .000001   .000001\nIMPROVE\nCALL FCN        3.\nEXIT\n")
            f.truncate()



