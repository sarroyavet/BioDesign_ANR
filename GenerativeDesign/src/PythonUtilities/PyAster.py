# Libraries {{{
import os
import sys
import numpy as np
from scipy.integrate import simps as simps
from scipy.spatial import cKDTree

# From code aster
from code_aster.Cata.Syntax import *
from code_aster.Cata.DataStructure import *
from code_aster.Cata.Commons import *
from Utilitai.partition import *

# In-house modules
sys.path.append('#srcDir/PythonUtilities')
from FemMesh2D import FemMesh2D as fm2
from FemVtk import FemVtk as fvt
import MorphoDesignFunctions as mdf
import meshMaker as mm
# }}}

# Macros {{{
# Read mesh {{{
def read_mesh(self, MESH_NAME, FORMAT, UNITE):
    ier = 0   # You need to return an integer
    # Add 1 for enumeration of the commands as the macro counts for one.
    self.set_icmd(1)

   # Declare "mesh" as the main output of the macro.
    self.DeclareOut("mesh", self.sd)

    # Import the operators/macros/procedures to use in your MACRO
    LIRE_MAILLAGE = self.get_cmd("LIRE_MAILLAGE")
    DEFI_FICHIER = self.get_cmd("DEFI_FICHIER")
    DEFI_FICHIER(ACTION = 'ASSOCIER',
            FICHIER = MESH_NAME,
            UNITE = UNITE)
    # This is the output of your macro
    mesh = LIRE_MAILLAGE(FORMAT = FORMAT,
                         UNITE = UNITE)
    DEFI_FICHIER(ACTION = 'LIBERER', UNITE = UNITE)
    return ier


READ_MESH = MACRO(
        nom="READ_MESH",
        op=read_mesh,
        sd_prod=maillage_sdaster,
        docu="Macro command that takes a mesh name as an input and returns a mesh object/concept",
        reentrant="n",
        fr="",
        MESH_NAME=SIMP(statut="o", typ='TXM'),
        FORMAT=SIMP(statut="d", typ='TXM', defaut = 'MED'),
        UNITE=SIMP(statut='d', typ='I', defaut=101),
)
# }}}

# Save results rmed {{{
def save_results_rmed(self, FILE_NAME, RESU, UNITE):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    # Add required code_aster functions
    DEFI_FICHIER = self.get_cmd("DEFI_FICHIER")
    IMPR_RESU = self.get_cmd("IMPR_RESU")
    _F = self.get_cmd("_F")
    # }}}
    # Macro body {{{
    DEFI_FICHIER(ACTION = 'ASSOCIER',
            FICHIER = FILE_NAME,
            UNITE = UNITE)
    IMPR_RESU(UNITE = UNITE,
              FORMAT = 'MED',
              RESU = RESU)
    DEFI_FICHIER(ACTION = 'LIBERER', UNITE = UNITE)
    # }}}
    return ier
# Macro heading {{{
SAVE_RESULTS_RMED = MACRO(
        nom = "SAVE_RESULTS_RMED",
        docu = "Macro to save results in rmed format indicating the output file name.",
        op = save_results_rmed,
        FILE_NAME = SIMP(statut = "o", typ = 'TXM'),
        RESU = FACT(statut='o', max='**',
                    regles = (AU_MOINS_UN('CHAM_GD', 'RESULTAT',
                                          'MAILLAGE'),
                              EXCLUS('CHAM_GD', 'RESULTAT'),
                              EXCLUS('TOUT_CMP', 'NOM_CMP')),
                    MAILLAGE = SIMP(statut = 'f',
                                    typ = (maillage_sdaster,
                                           squelette)),
                    CARA_ELEM = SIMP(statut = 'f',
                                     typ = cara_elem),
                    CHAM_GD = SIMP(statut ='f', typ = cham_gd_sdaster),
                    RESULTAT  = SIMP(statut = 'f',
                                     typ = resultat_sdaster),
                    INFO_MAILLAGE = SIMP(statut = 'f', typ = 'TXM',
                                         defaut = "NON",
                                         into = ("OUI", "NON"))),
        UNITE = SIMP(statut = 'd', typ = 'I', defaut = 101),
)
# }}}
# }}}

# Read two--geometry Mortar mesh {{{
def two_geometry_mortar(self, NAME_MAS, NAME_ESL, FORMAT_MAS, FORMAT_ESL,
        UNITE, GROUP_MA_MAIT, GROUP_MA_ESCL):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut('mesh', self.sd)
    # Add required code_aster functions
    ASSE_MAILLAGE = self.get_cmd("ASSE_MAILLAGE")
    MODI_MAILLAGE = self.get_cmd("MODI_MAILLAGE")
    _F = self.get_cmd("_F")
    CREA_MAILLAGE = self.get_cmd("CREA_MAILLAGE")
    DEFI_GROUP = self.get_cmd("DEFI_GROUP")
    # }}}
    # Macro body {{{
    # Read and assemble meshes
    __mesh_M = READ_MESH(MESH_NAME = NAME_MAS, FORMAT = FORMAT_MAS, UNITE = UNITE)
    __mesh_E = READ_MESH(MESH_NAME = NAME_ESL, FORMAT = FORMAT_ESL, UNITE = UNITE)
    __mesh0 = ASSE_MAILLAGE(MAILLAGE_1 = __mesh_M,
                            MAILLAGE_2 = __mesh_E,
                            OPERATION = 'SUPERPOSE')
    # Reorient contact boundaries
    __mesh0 = MODI_MAILLAGE(reuse = __mesh0,
                            MAILLAGE = __mesh0,
                            ORIE_PEAU_2D = _F(GROUP_MA = (GROUP_MA_MAIT,
                                                          GROUP_MA_ESCL)))
    # Create Mortar mesh
    mesh = CREA_MAILLAGE(MAILLAGE = __mesh0,
                         DECOUPE_LAC = _F(GROUP_MA_ESCL = (GROUP_MA_ESCL)))
    # Create ordered nodal groups form master and slave element groups
    mesh = DEFI_GROUP(CREA_GROUP_NO = (_F(GROUP_MA = (GROUP_MA_ESCL),
                                          NOM = GROUP_MA_ESCL,
                                          OPTION = 'NOEUD_ORDO')),
                      MAILLAGE = mesh,
                      reuse = mesh)
    mesh = DEFI_GROUP(CREA_GROUP_NO = (_F(GROUP_MA = (GROUP_MA_MAIT),
                                          NOM = GROUP_MA_MAIT,
                                          OPTION = 'NOEUD_ORDO')),
                      MAILLAGE = mesh,
                      reuse = mesh)
    # }}}
    return ier
TWO_GEOMETRY_MORTAR = MACRO(
        nom = "TWO_GEOMETRY_MORTAR",
        op = two_geometry_mortar,
        sd_prod = maillage_sdaster,
        docu = "Macro that returns a Mortar mesh from to mesh file names.",
        reentrant = "n",
        fr = "",
        NAME_MAS = SIMP(statut="o", typ='TXM'),
        NAME_ESL = SIMP(statut="o", typ='TXM'),
        GROUP_MA_MAIT = SIMP(statut="o", typ='TXM'),
        GROUP_MA_ESCL = SIMP(statut="o", typ='TXM'),
        FORMAT_MAS = SIMP(statut="d", typ='TXM', defaut = 'MED'),
        FORMAT_ESL = SIMP(statut="d", typ='TXM', defaut = 'MED'),
        UNITE = SIMP(statut='d', typ='I', defaut=101),
        )
# }}}

# Solve contact analysis by force {{{
def solve_contact_by_force(self, MESH, DIR_X, DIR_Y, RELOC_H,
        RELOC_GROUP_MA, GROUP_MA_MAIT, GROUP_MA_ESCL, CHAM_MATER,
        maxNewton, CONTACT, EXCIT, LIST_INST, MODELE, PREVRESU,
        CARA_ELEM):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut('resu', self.sd)
    # Add required code_aster functions
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    MODI_MAILLAGE = self.get_cmd("MODI_MAILLAGE")
    _F = self.get_cmd("_F")
    STAT_NON_LINE = self.get_cmd("STAT_NON_LINE")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    # }}}
    # Macro body {{{
    # Get the directed sign distance between the slave and master curves {{{
    if abs(np.sqrt((DIR_X**2.0 + DIR_Y**2.0)) - 1.0) > 1.0e-3:
        raise("Error: the direction vector is not a unit vector.")
    __mail_py = MAIL_PY()
    __mail_py.FromAster(MESH)
    meshFEM2D = fm2.MeshFromMail_Py(__mail_py)
    distance = meshFEM2D.DistanceTwoElementSets(GROUP_MA_MAIT + '_lines',
            GROUP_MA_ESCL + '_lines', direction = np.array([DIR_X, DIR_Y]))
    # }}}
    # Force the initial interpenetration {{{
    dx = (distance + RELOC_H)*DIR_X
    dy = (distance + RELOC_H)*DIR_Y
    __reloc1 = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                       OPERATION = 'AFFE',
                       MAILLAGE = MESH,
                       AFFE = (_F(TOUT = 'OUI',
                                  NOM_CMP = ('DX', 'DY'),
                                  VALE = (0.0, 0.0)),
                               _F(GROUP_MA = RELOC_GROUP_MA,
                                  NOM_CMP = ('DX', 'DY'),
                                  VALE = (dx, dy))))
    MODI_MAILLAGE(reuse = MESH,
                  MAILLAGE = MESH,
                  DEFORME = (_F(OPTION = 'TRAN',
                                DEPL = __reloc1)))
    # }}}
    # Solve the nonlinear problem {{{
    if PREVRESU == None:
        resu = STAT_NON_LINE(CHAM_MATER = CHAM_MATER,
                             COMPORTEMENT=_F(ITER_INTE_MAXI = maxNewton),
                             CONTACT = CONTACT,
                             CONVERGENCE=_F(ITER_GLOB_MAXI = maxNewton),
                             EXCIT = EXCIT,
                             CARA_ELEM = CARA_ELEM,
                             INFO = 1,
                             INCREMENT = _F(LIST_INST = LIST_INST,
                                            PRECISION = 1.E-06),
                             METHODE='NEWTON',
                             MODELE = MODELE,
                             NEWTON=_F(PREDICTION = 'TANGENTE'),
                             ARCHIVAGE=_F(PRECISION=1.E-06,
                                          CRITERE='RELATIF',))
    else:
        resu = STAT_NON_LINE(CHAM_MATER = CHAM_MATER,
                             COMPORTEMENT=_F(ITER_INTE_MAXI = maxNewton),
                             CONTACT = CONTACT,
                             CONVERGENCE=_F(ITER_GLOB_MAXI = maxNewton),
                             EXCIT = EXCIT,
                             INFO = 1,
                             INCREMENT = _F(LIST_INST = LIST_INST,
                                            PRECISION = 1.E-06),
                             METHODE='NEWTON',
                             MODELE = MODELE,
                             CARA_ELEM = CARA_ELEM,
                             ETAT_INIT = _F(EVOL_NOLI = PREVRESU,
                                            INST_ETAT_INIT = 0.0),
                             NEWTON=_F(PREDICTION = 'TANGENTE'),
                             ARCHIVAGE=_F(PRECISION=1.E-06,
                                          CRITERE='RELATIF',))
    # }}}
    ## Relocate the mesh {{{
    #__reloc2 = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
    #                      OPERATION = 'AFFE',
    #                      MAILLAGE = MESH,
    #                      AFFE = (_F(TOUT = 'OUI',
    #                                 NOM_CMP = ('DX', 'DY'),
    #                                 VALE = (0.0, 0.0)),
    #                           _F(GROUP_MA = RELOC_GROUP_MA,
    #                              NOM_CMP = ('DX', 'DY'),
    #                              VALE = (-dx, -dy))))
    #MODI_MAILLAGE(reuse = MESH,
    #              MAILLAGE = MESH,
    #              DEFORME = (_F(OPTION = 'TRAN',
    #                            DEPL = __reloc2)))
    ## }}}
    # Generate stress fields {{{
    CALC_CHAMP(reuse = resu,
               CONTRAINTE = ('SIGM_NOEU', 'SIGM_ELGA', 'SIEF_NOEU'),
               CRITERES = ('SIEQ_NOEU', 'SIEQ_ELGA'),
               FORCE = 'REAC_NODA',
               RESULTAT = resu)
    # }}}
    # }}}
    return ier
# Macro heading {{{
SOLVE_CONTACT_BY_FORCE = MACRO(
        nom = "SOLVE_CONTACT_BY_FORCE",
        op = solve_contact_by_force,
        sd_prod = evol_noli,
        docu = "Macro that returns the solution of a contact problem by force adding initial mesh penetration.",
        reentrant = "n",
        fr = "",
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        DIR_X = SIMP(statut = "o", typ = "R"),
        DIR_Y = SIMP(statut = "o", typ = "R"),
        RELOC_H = SIMP(statut = "o", typ = "R"),
        RELOC_GROUP_MA = SIMP(statut="o", typ='TXM'),
        GROUP_MA_MAIT = SIMP(statut="o", typ='TXM'),
        GROUP_MA_ESCL = SIMP(statut="o", typ='TXM'),
        CHAM_MATER = SIMP(statut = 'o', typ=cham_mater),
        maxNewton = SIMP(statut = "d", typ = "I", defaut = 10),
        CONTACT = SIMP(statut='o', typ=char_contact),
        EXCIT = FACT(statut = 'o', max='**',
                     CHARGE = SIMP(statut = 'o', typ = (char_meca,
                         char_cine_meca)),
                     FONC_MULT = SIMP(statut = 'f', typ = (fonction_sdaster,
                         nappe_sdaster, formule)),
                     TYPE_CHARGE = SIMP(statut='f', typ='TXM', defaut = "FIXE_CSTE",
                         into=("FIXE_CSTE", "FIXE_PILO", "SUIV",
                             "SUIV_PILO", "DIDI"))),
        PREVRESU  = SIMP(statut = "f", typ = evol_noli),
        LIST_INST = SIMP(statut = 'o', typ=(list_inst)),
        CARA_ELEM =SIMP(statut='f',typ=cara_elem),
        MODELE = SIMP(statut='o', typ=modele_sdaster)
        )
# }}}
# }}}

# Contact pressure as a function of ABSC_CURV {{{
def contact_pressure_absc_curv(self, RESU, MESH, GROUP_MA, GROUP_NO,
        INST, MODELE):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("pc", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CALC_PRESSION = self.get_cmd("CALC_PRESSION")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    POST_RELEVE_T = self.get_cmd("POST_RELEVE_T")
    RECU_FONCTION = self.get_cmd("RECU_FONCTION")
    CREA_TABLE = self.get_cmd("CREA_TABLE")
    # }}}
    # Macro body {{{
    # Set group_no
    if GROUP_NO == "__######__":
        group_no = GROUP_MA
    else:
        group_no = GROUP_NO
    # Extract the contact pressure as a function of ABSC_CURV
    if MODELE == None:
        __contPr = CALC_PRESSION(GEOMETRIE = 'INITIALE',
                                 GROUP_MA = (GROUP_MA),
                                 INST = INST,
                                 MAILLAGE = MESH,
                                 RESULTAT = RESU)
        strPression = 'LAGS_C'
    else:
        __contPr = CALC_PRESSION_FROM_FORMULA(
                MODELE = MODELE,
                RESU = RESU,
                NON_DEFORMED_MESH = MESH,
                GROUP_MA = GROUP_MA,
                INST = INST
                )
        strPression = 'X1'
    __PostR = POST_RELEVE_T(ACTION = (_F(INTITULE = 'Contact ' +
                                                    'pressure',
                                         OPERATION = 'EXTRACTION',
                                         GROUP_NO = group_no,
                                         CHAM_GD = __contPr,
                                         NOM_CMP = (strPression))))
    __contFo = RECU_FONCTION(PARA_X = 'ABSC_CURV',
                             PARA_Y = strPression,
                             TABLE = __PostR)
    pc = CREA_TABLE(FONCTION = _F(FONCTION = __contFo))
    # }}}
    return ier
# Macro heading {{{
CONTACT_PRESSURE_ABSC_CURV = MACRO(
        nom = 'CONTACT_PRESSURE_ABSC_CURV',
        op = contact_pressure_absc_curv,
        sd_prod = table_sdaster,
        docu = "Macro to compute the contact pressure as a function of ABSC_CURV.",
        MODELE = SIMP(statut='f', typ=modele_sdaster),
        RESU = SIMP(statut = "o", typ = evol_noli),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        GROUP_MA = SIMP(statut = "o", typ = 'TXM'),
        GROUP_NO = SIMP(statut = "d", typ = 'TXM', defaut = "__######__"),
        INST = SIMP(statut = "o", typ = "R")
)
# }}}
# }}}

# Compute mesh wear displacement {{{
def compute_mesh_wear_displacement(self, MODELE, GROUP_MA_CONT, RESU,
        MESH, INST, KAPPA, GROUP_NO_DX_0, GROUP_NO_DY_0, CHAM_MATER,
        CARA_ELEM):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("disp", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    AFFE_CHAR_MECA = self.get_cmd("AFFE_CHAR_MECA")
    MECA_STATIQUE = self.get_cmd("MECA_STATIQUE")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    # }}}
    # Macro body {{{
    __normal = CREA_CHAMP(TYPE_CHAM = 'NOEU_GEOM_R',
                          OPERATION = 'NORMALE',
                          MODELE = MODELE,
                          GROUP_MA = (GROUP_MA_CONT))
    try:
        __tbcont = CONTACT_PRESSURE_ABSC_CURV(RESU = RESU,
                                         MESH = MESH,
                                         GROUP_MA = GROUP_MA_CONT,
                                         INST = INST)
        # Get pressure vector
        p_c = np.array(__tbcont.EXTR_TABLE().values()['LAGS_C'])
    except:
        __tbcont = CONTACT_PRESSURE_ABSC_CURV(RESU = RESU,
                                         MESH = MESH,
                                         GROUP_MA = GROUP_MA_CONT,
                                         MODELE = MODELE,
                                         INST = INST)
        # Get pressure vector
        p_c = np.array(__tbcont.EXTR_TABLE().values()['X1'])
    # Compute the wear status
    __mail_py = MAIL_PY()
    __mail_py.FromAster(MESH)
    wear_status = mdf.WearDepth(p_c, __mail_py, GROUP_MA_CONT, __normal,
            KAPPA)
    # Set wear displacements as boundary conditions
    FFF1 = NodalWearDisplacements_F(wear_status, __mail_py)
    FFF2 = [_F(DX = 0.0, GROUP_NO = GROUP_NO_DX_0),
            _F(DY = 0.0, GROUP_NO = GROUP_NO_DY_0)]
    FFF = tuple(FFF1 + FFF2)
    __wearF =  AFFE_CHAR_MECA(DDL_IMPO = (FFF),
                             MODELE = MODELE)
    __wearR = MECA_STATIQUE(CHAM_MATER = CHAM_MATER,
                           EXCIT = (_F(CHARGE = __wearF)),
                           CARA_ELEM = CARA_ELEM,
                           MODELE = MODELE)
    __wearR = CALC_CHAMP(reuse = __wearR,
                        CONTRAINTE = ('SIGM_NOEU', 'SIGM_ELGA', 'SIEF_NOEU'),
                        CRITERES = ('SIEQ_ELGA'),
                        RESULTAT = __wearR)
    disp = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                      OPERATION = 'EXTR',
                      RESULTAT = __wearR,
                      NOM_CHAM = 'DEPL')
    # }}}
    return ier
# Macro heading {{{
COMPUTE_MESH_WEAR_DISPLACEMENT = MACRO(
        nom = 'COMPUTE_MESH_WEAR_DISPLACEMENT',
        op = compute_mesh_wear_displacement,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute a displacement field trigger by wear.",
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        GROUP_MA_CONT = SIMP(statut = "o", typ = 'TXM'),
        GROUP_NO_DX_0 = SIMP(statut='f',typ=grma ,validators=NoRepeat(),max='**'),
        GROUP_NO_DY_0 = SIMP(statut='f',typ=grma ,validators=NoRepeat(),max='**'),
        KAPPA = SIMP(statut = "o", typ = "R"),
        INST = SIMP(statut = "d", typ = "R", defaut = 1.0),
        CHAM_MATER = SIMP(statut = 'o', typ=cham_mater),
        RESU = SIMP(statut = "o", typ = evol_noli),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        CARA_ELEM =SIMP(statut='f',typ=cara_elem),
        )
# }}}
# }}}

# Save results vtk {{{
def save_results_vtk(self, FILE_NAME, RESU, MESH, GROUP_MA_SURF, CHAMPS,
        INST, NUME_RESU_GROWTH, DIMS, NODE_FIELDS):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Add required code_aster functions
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    FORMULE = self.get_cmd("FORMULE")
    _F = self.get_cmd("_F")
    # }}}
    # Macro body {{{
    if GROUP_MA_SURF == None:
        group = []
    else:
        group = [GROUP_MA_SURF]
    # Necessary formulas and fields {{{
    __Shyd_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ'),
                     VALE = 'Shyd(SIXX, SIYY, SIZZ)',
                     Shyd = mdf.Shyd_pln_strain)
    __Svon_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ', 'SIXY'),
                     VALE = 'Svon(SIXX, SIYY, SIZZ, SIXY)',
                     Svon = mdf.Svmis_pln_strain)
    CALC_CHAMP(reuse = RESU,
               CONTRAINTE = ('SIGM_NOEU', 'SIGM_ELGA', 'SIEF_NOEU'),
               CRITERES = ('SIEQ_ELGA'),
               RESULTAT = RESU)
    # }}}
    __mail_py = MAIL_PY()
    __mail_py.FromAster(MESH)
    mesh_2d = fm2.MeshFromMail_Py(__mail_py, *group)
    vtkResults = fvt(mesh_2d)
    # Add fields to vtk {{{
    for champ_name in CHAMPS:
        if (champ_name == 'DEPL') or (champ_name == 'disp'):
        # Add displacement {{{
            __disp = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                                OPERATION = 'EXTR',
                                RESULTAT = RESU,
                                INST = INST,
                                NOM_CHAM = 'DEPL')
            disp_x = __disp.EXTR_COMP('DX', []).valeurs
            disp_y = __disp.EXTR_COMP('DY', []).valeurs
            if DIMS == 2:
                disp_z = np.zeros(disp_x.size)
            elif DIMS == 3:
                disp_z = __disp.EXTR_COMP('DZ', []).valeurs
            array = np.column_stack((disp_x, disp_y, disp_z))
        # }}}
        elif champ_name == 'STRESS2D':
            # Add stress 2D {{{
            __stress = CREA_CHAMP(TYPE_CHAM = 'NOEU_SIEF_R',
                                  OPERATION = 'EXTR',
                                  RESULTAT = RESU,
                                  INST = INST,
                                  NOM_CHAM = 'SIEF_NOEU')
            if DIMS == 2:
                stresses_tensor_xx = __stress.EXTR_COMP('SIXX', []).valeurs
                stresses_tensor_yy = __stress.EXTR_COMP('SIYY', []).valeurs
                stresses_tensor_zz = __stress.EXTR_COMP('SIZZ', []).valeurs
                stresses_tensor_xy = __stress.EXTR_COMP('SIXY', []).valeurs
                array = np.zeros((len(stresses_tensor_xy), 3, 3))
                for k1 in range(len(stresses_tensor_xy)):
                    array[k1, 0, 0] = stresses_tensor_xx[k1]
                    array[k1, 1, 1] = stresses_tensor_yy[k1]
                    array[k1, 0, 1] = stresses_tensor_xy[k1]
                    array[k1, 1, 0] = stresses_tensor_xy[k1]
                    array[k1, 2, 2] = stresses_tensor_zz[k1]
            else:
                stresses_tensor = __stress.EXTR_COMP('', []).valeurs
                array = fvt.FromVoigtToTensor(stresses_tensor, '2D')
            # }}}
        elif (champ_name == 'VONMISES') or (champ_name == 'tauMis'):
            # Add von Mises stress {{{
            CALC_CHAMP(reuse = RESU,
                       CONTRAINTE = ('SIGM_NOEU', 'SIGM_ELGA',
                                     'SIEF_NOEU'),
                       CRITERES = ('SIEQ_ELGA'),
                       INST = INST,
                       RESULTAT = RESU)
            __resu = CALC_CHAMP(CHAM_UTIL = _F(FORMULE = (__Svon_f),
                                NOM_CHAM = 'SIGM_NOEU',
                                NUME_CHAM_RESU = 2),
                                INST = INST,
                                RESULTAT = RESU)
            __vmis = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                                OPERATION = 'EXTR',
                                RESULTAT = __resu,
                                INST = INST,
                                NOM_CHAM = 'UT02_NOEU')
            array = __vmis.EXTR_COMP('X1', []).valeurs
            # }}}
        elif (champ_name == 'SHYD') or (champ_name == 'sigHyd'):
            # Add hydrostatic stress {{{
            CALC_CHAMP(reuse = RESU,
                       CONTRAINTE = ('SIGM_NOEU', 'SIGM_ELGA',
                                     'SIEF_NOEU'),
                       CRITERES = ('SIEQ_ELGA'),
                       INST = INST,
                       RESULTAT = RESU)
            __resu = CALC_CHAMP(CHAM_UTIL = _F(FORMULE = (__Shyd_f),
                                NOM_CHAM = 'SIGM_NOEU',
                                NUME_CHAM_RESU = 2),
                                INST = INST,
                                RESULTAT = RESU)
            __shyd = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                                OPERATION = 'EXTR',
                                RESULTAT = __resu,
                                INST = INST,
                                NOM_CHAM = 'UT02_NOEU')
            array = __shyd.EXTR_COMP('X1', []).valeurs
            # }}}
        elif champ_name == 'growthData':
            # Add growth function {{{
            nomcham = 'UT' + str(NUME_RESU_GROWTH).zfill(2) + '_NOEU'
            __grw_fi = CREA_CHAMP(OPERATION = 'EXTR',
                                  TYPE_CHAM = 'NOEU_NEUT_R',
                                  NOM_CHAM = nomcham,
                                  INST = INST,
                                  RESULTAT = RESU)
            grwStr = __grw_fi.EXTR_COMP('', [], topo = 1)
            array = np.zeros(mesh_2d.nodes.numNods)
            for k1 in range(grwStr.valeurs.size):
                array[grwStr.noeud[k1] - 1] = grwStr.valeurs[k1]
            # }}}
        else:
            # Otherwise error {{{
            raise("Error: unknown keyword in CHAMPS.")
            # }}}
        vtkResults.add_data(array, champ_name)
    # }}}
    # Add field to vtk from NODE_FIELDS {{{
    if NODE_FIELDS != None:
        for node_field_i in NODE_FIELDS.List_F():
            champ_name = node_field_i["NOM"]
            nom_cmp_i = node_field_i["NOM_CMP"]
            cham_i = node_field_i["CHAM_NO"]
            array = cham_i.EXTR_COMP(nom_cmp_i, []).valeurs
            vtkResults.add_data(array, champ_name)
    # }}}
    vtkResults.WriteVtk(FILE_NAME)
    # }}}
    return ier
# Macro heading {{{
SAVE_RESULTS_VTK = MACRO(
        nom = "SAVE_RESULTS_VTK",
        docu = "Macro to save results in vtk format indicating the output file name.",
        op = save_results_vtk,
        FILE_NAME = SIMP(statut = "o", typ = 'TXM'),
        RESU = SIMP(statut = "o", typ = resultat_sdaster),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        GROUP_MA_SURF = SIMP(statut="f", typ='TXM'),
        CHAMPS = SIMP(statut='f',typ=grma ,validators=NoRepeat(),max='**'),
        INST = SIMP(statut = "f", typ = "R"),
        DIMS = SIMP(statut = "d", typ = "I", defaut = 2, into=( 2 , 3)),
        NUME_RESU_GROWTH = SIMP(statut = "d", typ = "I", defaut = 10, val_min = 2),
        NODE_FIELDS = FACT(statut = 'f', max = '**',
                           CHAM_NO = SIMP(statut = 'o',
                                          typ = cham_no_sdaster),
                           NOM = SIMP(statut = "o", typ = "TXM"),
                           NOM_CMP = SIMP(statut = "o", typ = "TXM"))
)
# }}}
# }}}

# Morphogenesis growth stress {{{
def morphogenesis_growth(self, RESU, INST, NUME_RESU_GROWTH, GROUP_MA,
        ALPHA, SIGLIM, TAULIM, VEL, MODELE, MIN_FUN, MAX_FUN, XMIN,
        XMAX, YMIN, YMAX, MESH, GFUNC, COEFF):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("growth_t", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    FORMULE = self.get_cmd("FORMULE")
    # }}}
    # Macro body {{{
    # Test input {{{
    if MIN_FUN >= MAX_FUN:
        raise("Error: MIN_FUN must be lower than MAX_FUN")
    # }}}
    # Compute necessary stress fields {{{
    CALC_CHAMP(reuse = RESU,
               CONTRAINTE = ('SIGM_NOEU', 'SIGM_ELGA',
                             'SIEF_NOEU'),
               CRITERES = ('SIEQ_ELGA'),
               INST = INST,
               RESULTAT = RESU)
    # }}}
    # Compute hydrostatic and shear stress {{{
    # Formulas
    __Shyd_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ'),
                     VALE = 'Shyd(SIXX, SIYY, SIZZ)',
                     Shyd = mdf.Shyd_pln_strain)
    __Svon_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ', 'SIXY'),
                     VALE = 'Svon(SIXX, SIYY, SIZZ, SIXY)',
                     Svon = mdf.Svmis_pln_strain)
    # Stresses as result
    nume1 = NUME_RESU_GROWTH - 1
    CALC_CHAMP(CHAM_UTIL = _F(FORMULE = (__Shyd_f, __Svon_f),
                              NOM_CHAM = 'SIGM_NOEU',
                              NUME_CHAM_RESU = nume1),
               RESULTAT = RESU,
               INST = INST,
               reuse = RESU)
    # Stresses as fields
    __vonStr = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                          OPERATION = 'EXTR',
                          RESULTAT = RESU,
                          INST = INST,
                          NOM_CHAM = 'UT' + str(nume1).zfill(2) + '_NOEU')
    # Box limits {{{
    # To select the region from which the maximum values will be
    # computed
    # BoxIn function {{{
    def BoxIn(X, Y, argDict):
        if "xmin" in argDict:
            byXmin = X >= argDict["xmin"]
        else:
            byXmin = True
        if "xmax" in argDict:
            byXmax = X <= argDict["xmax"]
        else:
            byXmax = True
        if "ymin" in argDict:
            byYmin = Y >= argDict["ymin"]
        else:
            byYmin = True
        if "ymax" in argDict:
            byYmax = Y <= argDict["ymax"]
        else:
            byYmax = True
        byX = byXmin and byYmax
        byY = byYmin and byYmax
        if byX and byY:
            return 1.0
        else:
            return 0.0
    # }}}
    # Make limits {{{
    limits = {}
    if XMIN:
        limits["xmin"] = XMIN
    if XMAX:
        limits["xmax"] = XMAX
    if YMIN:
        limits["ymin"] = YMIN
    if YMAX:
        limits["ymax"] = YMAX
    # }}}
    __geom = CREA_CHAMP(OPERATION = 'EXTR',
                        TYPE_CHAM = 'NOEU_GEOM_R',
                        NOM_CHAM = 'GEOMETRIE',
                        MAILLAGE = MESH)
    __boxfo1 = FORMULE(NOM_PARA = ("X", "Y", "X1"),
                       VALE = ('X1*func(X, Y, argDict)'),
                       func = BoxIn,
                       argDict = limits)
    __boxfo2 = FORMULE(NOM_PARA = ("X", "Y", "X2"),
                       VALE = ('X2*func(X, Y, argDict)'),
                       func = BoxIn,
                       argDict = limits)
    __box_fu = CREA_CHAMP(OPERATION = 'AFFE',
            TYPE_CHAM = 'NOEU_NEUT_F',
            MODELE = MODELE,
            AFFE = _F(TOUT = 'OUI',
                      NOM_CMP = ('X1', 'X2'),
                      VALE_F = (__boxfo1, __boxfo2)))
    __box_ev = CREA_CHAMP(OPERATION = 'EVAL',
                          TYPE_CHAM = 'NOEU_NEUT_R',
                          CHAM_F = __box_fu,
                          CHAM_PARA = (__geom, __vonStr))
    # }}}
    # }}}
    # Compute growth stress tensor {{{
    # Maximum stresses
    tauMiswp = __box_ev.EXTR_COMP('X2', [GROUP_MA]).valeurs
    sigHydwp = __box_ev.EXTR_COMP('X1', [GROUP_MA]).valeurs
    maxSigHydCom = abs(min(sigHydwp))
    maxTau = max(tauMiswp)
    # Formula
    sgfunc = eval(GFUNC)
    if COEFF == None:
        func = lambda X1, X2: sgfunc(X1/maxSigHydCom, X2/maxTau, ALPHA,
                shrlim = TAULIM, hydlim = SIGLIM, vel = VEL,
                minFun = MIN_FUN, maxFun = MAX_FUN)
    else:
        func = lambda X1, X2: sgfunc(X1/maxSigHydCom, X2/maxTau, ALPHA,
                COEFF,
                shrlim = TAULIM, hydlim = SIGLIM, vel = VEL,
                minFun = MIN_FUN, maxFun = MAX_FUN)
    __Sgrw_f  = FORMULE(NOM_PARA = ('X1', 'X2'),
                        VALE = 'func(X1, X2)',
                        func = func,
                        )
    #  Growth stress {as result}
    nume2 = NUME_RESU_GROWTH
    CALC_CHAMP(CHAM_UTIL = (_F(FORMULE = (__Sgrw_f),
                              NOM_CHAM = 'UT' + str(nume1).zfill(2) + '_NOEU',
                              NUME_CHAM_RESU = nume2)),
               RESULTAT = RESU,
               INST = INST,
               reuse = RESU)
    #  Growth stress {as field}
    __grw_fi = CREA_CHAMP(OPERATION = 'EXTR',
                          TYPE_CHAM = 'NOEU_NEUT_R',
                          NOM_CHAM = 'UT' + str(nume2).zfill(2) + '_NOEU',
                          INST = INST,
                          RESULTAT = RESU)
    # Growth stress as tensor
    growth_t = CREA_CHAMP(OPERATION = 'ASSE',
                          TYPE_CHAM = 'NOEU_SIEF_R',
                          MODELE = MODELE,
                          ASSE = (_F(GROUP_MA = GROUP_MA,
                                    CHAM_GD = __grw_fi,
                                    NOM_CMP = ('X1'),
                                    NOM_CMP_RESU = ('SIXX')),
                                  _F(GROUP_MA = GROUP_MA,
                                    CHAM_GD = __grw_fi,
                                    NOM_CMP = ('X1'),
                                    NOM_CMP_RESU = ('SIYY')),
                                  _F(GROUP_MA = GROUP_MA,
                                    CHAM_GD = __grw_fi,
                                    NOM_CMP = ('X1'),
                                    NOM_CMP_RESU = ('SIZZ')),),)
    # }}}
    # }}}
    return ier
# Macro heading {{{
MORPHOGENESIS_GROWTH = MACRO(
        nom = 'MORPHOGENESIS_GROWTH',
        op = morphogenesis_growth,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute the morphogenesis-inspired growth stress tensor.",
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        GROUP_MA = SIMP(statut = "o", typ = 'TXM'),
        ALPHA = SIMP(statut = "o", typ = "R"),
        SIGLIM = SIMP(statut = "o", typ = "R"),
        TAULIM = SIMP(statut = "o", typ = "R"),
        VEL = SIMP(statut = "o", typ = "R"),
        INST = SIMP(statut = "f", typ = "R"),
        MIN_FUN = SIMP(statut = "d", typ = "R", defaut = 0.0),
        MAX_FUN = SIMP(statut = "d", typ = "R", defaut = 1.0),
        XMIN = SIMP(statut = "f", typ = "R"),
        XMAX = SIMP(statut = "f", typ = "R"),
        YMIN = SIMP(statut = "f", typ = "R"),
        YMAX = SIMP(statut = "f", typ = "R"),
        NUME_RESU_GROWTH = SIMP(statut = "d", typ = "I", defaut = 10, val_min = 2),
        RESU = SIMP(statut = "o", typ = resultat_sdaster),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        GFUNC = SIMP(statut = "d", typ = "TXM", defaut = 'mdf.Sgrowth'),
        COEFF = SIMP(statut = "f", typ = "R"),
        )
# }}}
# }}}

# Compute mesh morphogenesis displacement {{{
def compute_mesh_morphogenesis_displacement(self, RESU, INST,
        NUME_RESU_GROWTH, GROUP_MA, ALPHA, SIGLIM, TAULIM, VEL, MODELE,
        EXCIT, CHAM_MATER, REF_FIXED_DX, REF_FIXED_DY, CARA_ELEM,
        MIN_FUN, MAX_FUN, MESH, XMIN, XMAX, YMIN, YMAX, GFUNC, COEFF):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("disp", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    AFFE_CHAR_MECA = self.get_cmd("AFFE_CHAR_MECA")
    MECA_STATIQUE = self.get_cmd("MECA_STATIQUE")
    STAT_NON_LINE = self.get_cmd("STAT_NON_LINE")
    FORMULE = self.get_cmd("FORMULE")
    # }}}
    # Macro body {{{
    # Compute growth tensor {{{
    __grw_t = MORPHOGENESIS_GROWTH(
            RESU = RESU,
            INST = INST,
            NUME_RESU_GROWTH = NUME_RESU_GROWTH,
            GROUP_MA = GROUP_MA,
            ALPHA = ALPHA,
            SIGLIM = SIGLIM,
            TAULIM = TAULIM,
            VEL = VEL,
            MODELE = MODELE,
            MIN_FUN = MIN_FUN,
            MAX_FUN = MAX_FUN,
            MESH = MESH,
            XMIN = XMIN,
            XMAX = XMAX,
            YMIN = YMIN,
            YMAX = YMAX,
            COEFF = COEFF,
            GFUNC = GFUNC,
            )
    # }}}
    # Compute displacement {{{
    disp = COMPUTE_MESH_DISPLACEMENT_FROM_GROWTH_FORCE(
            GRW_TEN = __grw_t,
            MODELE = MODELE,
            EXCIT = EXCIT,
            CHAM_MATER = CHAM_MATER,
            CARA_ELEM = CARA_ELEM,
            REF_FIXED_DX = REF_FIXED_DX,
            REF_FIXED_DY = REF_FIXED_DY,
            GROUP_MA = GROUP_MA)
    # }}}
    # }}}
    return ier
# Macro heading {{{
COMPUTE_MESH_MORPHOGENESIS_DISPLACEMENT = MACRO(
        nom = 'COMPUTE_MESH_MORPHOGENESIS_DISPLACEMENT',
        op = compute_mesh_morphogenesis_displacement,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute the morphogenesis-inspired displacement.",
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        GROUP_MA = SIMP(statut = "o", typ = 'TXM'),
        REF_FIXED_DX = SIMP(statut='f',typ=grma,validators=NoRepeat(),max='**', min = 0),
        REF_FIXED_DY = SIMP(statut='f',typ=grma,validators=NoRepeat(),max='**', min = 0),
        ALPHA = SIMP(statut = "o", typ = "R"),
        SIGLIM = SIMP(statut = "o", typ = "R"),
        TAULIM = SIMP(statut = "o", typ = "R"),
        VEL = SIMP(statut = "o", typ = "R"),
        INST = SIMP(statut = "f", typ = "R"),
        NUME_RESU_GROWTH = SIMP(statut = "d", typ = "I", defaut = 10),
        RESU = SIMP(statut = "o", typ = resultat_sdaster),
        CHAM_MATER = SIMP(statut = 'o', typ=cham_mater),
        EXCIT = FACT(statut = 'o', max='**',
                     CHARGE = SIMP(statut = 'o', typ = (char_meca,
                         char_cine_meca)),
                     FONC_MULT = SIMP(statut = 'f', typ = (fonction_sdaster,
                         nappe_sdaster, formule)),
                     TYPE_CHARGE = SIMP(statut='f', typ='TXM', defaut = "FIXE_CSTE",
                         into=("FIXE_CSTE", "FIXE_PILO", "SUIV",
                             "SUIV_PILO", "DIDI"))),
        CARA_ELEM =SIMP(statut='f',typ=cara_elem),
        MIN_FUN = SIMP(statut = "d", typ = "R", defaut = 0.0),
        MAX_FUN = SIMP(statut = "d", typ = "R", defaut = 1.0),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        XMIN = SIMP(statut = "f", typ = "R"),
        XMAX = SIMP(statut = "f", typ = "R"),
        YMIN = SIMP(statut = "f", typ = "R"),
        YMAX = SIMP(statut = "f", typ = "R"),
        GFUNC = SIMP(statut = "d", typ = "TXM", defaut = 'mdf.Sgrowth'),
        COEFF = SIMP(statut = "f", typ = "R"),
        )
# }}}
# }}}

# Compute mesh morphogenesis displacement with contact {{{
def compute_mesh_morphogenesis_displacement_with_contact(self, RESU, INST,
        NUME_RESU_GROWTH, GROUP_MA, ALPHA, SIGLIM, TAULIM, VEL, MODELE,
        EXCIT, CHAM_MATER, minSteps, maxSteps, maxNewton, CONTACT,
        SIMULTANEOUS, REF_FIXED_GROUP_NO, CARA_ELEM):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("disp", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    AFFE_CHAR_MECA = self.get_cmd("AFFE_CHAR_MECA")
    MECA_STATIQUE = self.get_cmd("MECA_STATIQUE")
    STAT_NON_LINE = self.get_cmd("STAT_NON_LINE")
    FORMULE = self.get_cmd("FORMULE")
    DEFI_FONCTION = self.get_cmd("DEFI_FONCTION")
    DEFI_LIST_REEL = self.get_cmd("DEFI_LIST_REEL")
    DEFI_LIST_INST = self.get_cmd("DEFI_LIST_INST")
    # }}}
    # Macro body {{{
    # Tests group inputs {{{
    if isinstance(GROUP_MA, str):
        groups = [GROUP_MA]
    else:
        groups = GROUP_MA
    if isinstance(REF_FIXED_GROUP_NO, str):
        fixedGroup = REF_FIXED_GROUP_NO
    else:
        if REF_FIXED_GROUP_NO == None:
            fixedGroup = None
        else:
            if len(REF_FIXED_GROUP_NO) == 1:
                fixedGroup = REF_FIXED_GROUP_NO[0]
            elif len(REF_FIXED_GROUP_NO) == 0:
                fixedGroup = None
            else:
                raise("Error: REF_FIXED_GROUP_NO includes several groups.")
    # }}}
    # Compute growth tensor {{{
    numGroups = len(groups)
    __gt = [None]*numGroups
    __gg = [None]*numGroups
    __st = [None]*numGroups
    for k1 in range(numGroups):
        __gt[k1] = MORPHOGENESIS_GROWTH(
                   RESU = RESU,
                   INST = INST,
                   NUME_RESU_GROWTH = NUME_RESU_GROWTH + k1*2,
                   GROUP_MA = groups[k1],
                   ALPHA = ALPHA,
                   SIGLIM = SIGLIM,
                   TAULIM = TAULIM,
                   VEL = VEL,
                   MODELE = MODELE
                   )
        __gg[k1] = CREA_CHAMP(OPERATION = 'DISC',
                              TYPE_CHAM = 'ELGA_SIEF_R',
                              MODELE = MODELE,
                              CHAM_GD = __gt[k1],
                              PROL_ZERO = 'OUI')
        __st[k1] = AFFE_CHAR_MECA(MODELE = MODELE,
                                  PRE_SIGM = _F(SIGM = __gg[k1]))
    # }}}
    # Make growth load and time list {{{
    loads = []
    for load in EXCIT:
        loads.append(load)
    __tm = [None]*numGroups
    t0 = INST
    intervalle_list = [_F(JUSQU_A = t0, NOMBRE = minSteps)]
    for k1 in range(numGroups):
        if SIMULTANEOUS == 1:
            vales = (t0, 0.0, t0 + 1.0, 1.0, t0 + 2.0, 0.0)
        elif SIMULTANEOUS == 2:
            vales = (t0, 0.0, t0 + 1.0, 1.0)
        __tm[k1] = DEFI_FONCTION(NOM_PARA='INST',
                                 VALE= vales,
                                 PROL_GAUCHE = 'CONSTANT',
                                 PROL_DROITE = 'CONSTANT')
        loads.append(_F(CHARGE = __st[k1], FONC_MULT = __tm[k1]))
        t0 += 1
        intervalle_list.append(_F(JUSQU_A = t0, NOMBRE = minSteps))
    __list = DEFI_LIST_REEL(DEBUT = 0.0, INTERVALLE = intervalle_list)
    __time = DEFI_LIST_INST(DEFI_LIST = _F(LIST_INST = __list),
                            ECHEC = _F(SUBD_METHODE = 'MANUEL',
                                       SUBD_PAS = 2,
                                       SUBD_PAS_MINI = 1.0/maxSteps,
                                       SUBD_NIVEAU = maxSteps))
    # }}}
    # Solve growth linear problem {{{
    __resu = STAT_NON_LINE(
            ETAT_INIT = _F(EVOL_NOLI = RESU),
            CHAM_MATER = CHAM_MATER,
            COMPORTEMENT=_F(ITER_INTE_MAXI = maxNewton),
            CONTACT = CONTACT,
            CONVERGENCE=_F(ITER_GLOB_MAXI = maxNewton),
            EXCIT = loads,
            INFO = 1,
            INCREMENT = _F(LIST_INST = __time,
                           PRECISION = 1.E-06),
            METHODE='NEWTON',
            MODELE = MODELE,
            CARA_ELEM = CARA_ELEM,
            NEWTON=_F(PREDICTION = 'TANGENTE'),
            ARCHIVAGE=_F(PRECISION=1.E-06,
                         CRITERE='RELATIF',)
            )
    # }}}
    # Compute growth displacement {{{
    # Formulas
    __dispxf = FORMULE(NOM_PARA = ('DX'), VALE = 'DX')
    __dispyf = FORMULE(NOM_PARA = ('DY'), VALE = 'DY')
    __dispze = FORMULE(NOM_PARA = ('DX','DY'), VALE = '0')
    __df = [None]*numGroups
    for k1 in range(numGroups):
        __df[k1] = CREA_CHAMP(OPERATION = 'AFFE',
                              TYPE_CHAM = 'NOEU_NEUT_F',
                              MODELE = MODELE,
                              AFFE = (_F(TOUT = 'OUI',
                                         NOM_CMP = ('X1', 'X2'),
                                         VALE_F = (__dispze, __dispze)),
                                      _F(GROUP_MA = groups[k1],
                                         NOM_CMP = ('X1', 'X2'),
                                         VALE_F = (__dispxf, __dispyf))))
    # Extraction of displacement for the reference instance INST, d_0
    __disp0 = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                         OPERATION = 'EXTR',
                         RESULTAT = __resu,
                         INST = INST,
                         NOM_CHAM = 'DEPL')
    if SIMULTANEOUS == 1:
        __dex = [None]*numGroups
        __dls = [None]*numGroups
        __dev = [None]*numGroups
        __dch = [None]*numGroups
        t0 = INST
        for k1 in range(numGroups):
            t0 += 1.0
            # Extraction of displacement for each instant d_i
            __dex[k1] = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                                   OPERATION = 'EXTR',
                                   RESULTAT = __resu,
                                   INST = t0,
                                   NOM_CHAM = 'DEPL')
            # Keep only the displacement produce by the new loads: d_i - d_0
            __dls[k1] = CREA_CHAMP(OPERATION = 'COMB',
                                   TYPE_CHAM = 'NOEU_DEPL_R',
                                   COMB = (_F(CHAM_GD = __dex[k1], COEF_R =  1.0),
                                           _F(CHAM_GD = __disp0, COEF_R = -1.0)))
            # Evaluate the __df: make zero __dls outside the specified group
            __dev[k1] = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                                   OPERATION = 'EVAL',
                                   CHAM_F = __df[k1],
                                   CHAM_PARA = __dls[k1])
            # Convert the evaluated (neutral) fields into displacement fields
            __dch[k1] = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                                   OPERATION = 'ASSE',
                                   MODELE = MODELE,
                                   ASSE = _F(CHAM_GD = __dev[k1],
                                             TOUT = 'OUI',
                                             NOM_CMP = ('X1', 'X2'),
                                             NOM_CMP_RESU = ('DX', 'DY')))
    elif SIMULTANEOUS == 2:
        # Extraction of displacement at the last instant: d_i
        __dex = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                           OPERATION = 'EXTR',
                           RESULTAT = __resu,
                           INST = INST + numGroups,
                           NOM_CHAM = 'DEPL')
        # Keep only the displacement produce by the new loads: d_i - d_0
        __dls = CREA_CHAMP(OPERATION = 'COMB',
                           TYPE_CHAM = 'NOEU_DEPL_R',
                           COMB = (_F(CHAM_GD = __dex, COEF_R =  1.0),
                                   _F(CHAM_GD = __disp0, COEF_R = -1.0)))
        __dev = [None]*numGroups
        __dch = [None]*numGroups
        for k1 in range(numGroups):
            # Evaluate the __df: make zero __dls outside the specified group
            __dev[k1] = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                                   OPERATION = 'EVAL',
                                   CHAM_F = __df[k1],
                                   CHAM_PARA = __dls)
            # Convert the evaluated (neutral) fields into displacement fields
            __dch[k1] = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                                   OPERATION = 'ASSE',
                                   MODELE = MODELE,
                                   ASSE = _F(CHAM_GD = __dev[k1],
                                             TOUT = 'OUI',
                                             NOM_CMP = ('X1', 'X2'),
                                             NOM_CMP_RESU = ('DX', 'DY')))
    # Sum up all the displacement fields
    combList = []
    for k1 in range(numGroups):
        combList.append(_F(CHAM_GD = __dch[k1], COEF_R =  1.0))
    if fixedGroup == None:
        # Add all the computed displacements and return
        disp = CREA_CHAMP(OPERATION = 'COMB',
                          TYPE_CHAM = 'NOEU_DEPL_R',
                          COMB = combList)
    else:
        # Add all the computed displacements
        __disp = CREA_CHAMP(OPERATION = 'COMB',
                            TYPE_CHAM = 'NOEU_DEPL_R',
                            COMB = combList)
        # Get the displacement assign to the reference node, the node
        # where the displacement should be zero
        refDx = __disp.EXTR_COMP('DX', [fixedGroup]).valeurs
        refDy = __disp.EXTR_COMP('DY', [fixedGroup]).valeurs
        if len(refDx) == 1:
            refDx = refDx[0]
            refDy = refDy[0]
        else:
            raise("Error: the reference group includes several nodes.")
        # Create a CONSTANT displacement field equal to minus the
        # Displacement of the reference field.
        __dispRf = CREA_CHAMP(OPERATION = 'AFFE',
                              TYPE_CHAM = 'NOEU_DEPL_R',
                              MODELE = MODELE,
                              AFFE = _F(GROUP_MA = groups,
                                        NOM_CMP = ('DX', 'DY'),
                                        VALE = (refDx, refDy)))
        # Compute __disp - __dispRf to obtain a displacement field where
        # the displacements at the reference point are zero.
        disp = CREA_CHAMP(OPERATION = 'COMB',
                          TYPE_CHAM = 'NOEU_DEPL_R',
                          COMB = (_F(CHAM_GD = __disp, COEF_R =  1.0),
                                  _F(CHAM_GD = __dispRf, COEF_R = -1.0)))
    # }}}
    # }}}
    return ier
# Macro heading {{{
COMPUTE_MESH_MORPHOGENESIS_DISPLACEMENT_WITH_CONTACT = MACRO(
        nom = 'COMPUTE_MESH_MORPHOGENESIS_DISPLACEMENT_WITH_CONTACT',
        op = compute_mesh_morphogenesis_displacement_with_contact,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute the morphogenesis-inspired displacement.",
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        GROUP_MA = SIMP(statut='o',typ=grma,validators=NoRepeat(),max='**'),
        REF_FIXED_GROUP_NO = SIMP(statut='f',typ=grma,validators=NoRepeat(),max='**', min = 0),
        ALPHA = SIMP(statut = "o", typ = "R"),
        SIGLIM = SIMP(statut = "o", typ = "R"),
        TAULIM = SIMP(statut = "o", typ = "R"),
        VEL = SIMP(statut = "o", typ = "R"),
        INST = SIMP(statut = "d", typ = "R", defaut = -1.0),
        NUME_RESU_GROWTH = SIMP(statut = "d", typ = "I", defaut = 10),
        RESU = SIMP(statut = "o", typ = evol_noli),
        CHAM_MATER = SIMP(statut = 'o', typ=cham_mater),
        EXCIT = FACT(statut = 'o', max='**',
                     CHARGE = SIMP(statut = 'o', typ = (char_meca,
                         char_cine_meca)),
                     FONC_MULT = SIMP(statut = 'f', typ = (fonction_sdaster,
                         nappe_sdaster, formule)),
                     TYPE_CHARGE = SIMP(statut='f', typ='TXM', defaut = "FIXE_CSTE",
                         into=("FIXE_CSTE", "FIXE_PILO", "SUIV",
                             "SUIV_PILO", "DIDI"))),
        maxNewton = SIMP(statut = "d", typ = "I", defaut = 10),
        minSteps = SIMP(statut = "d", typ = "I", defaut = 1),
        maxSteps = SIMP(statut = "d", typ = "I", defaut = 2),
        CONTACT = SIMP(statut='o', typ=char_contact),
        SIMULTANEOUS = SIMP(statut='d', typ='I', defaut= 1, into=( 1 , 2)),
        CARA_ELEM =SIMP(statut='f',typ=cara_elem),
        )
# }}}
# }}}

# Contact area from resultat cont_elem {{{
def contact_area_cont_elem(self, RESU, MESH, INST, MAX_JEU):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("table", self.sd)
    # Add required code_aster functions
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    CREA_TABLE = self.get_cmd("CREA_TABLE")
    _F = self.get_cmd("_F")
    # }}}
    # Macro body {{{
    # Make mesh 2d
    __mail_py = MAIL_PY()
    __mail_py.FromAster(MESH)
    # Extract contact results
    __contE = CREA_CHAMP(TYPE_CHAM = 'ELEM_CLAC_R',
                         OPERATION = 'EXTR',
                         RESULTAT = RESU,
                         INST = INST,
                         NOM_CHAM = 'CONT_ELEM')
    contE_t = __contE.EXTR_COMP('CONT', [], topo = 1)
    jeuE_t = __contE.EXTR_COMP('JEU', [], topo = 1)
    limits = mdf.Limits_active_CONT_ELEM1D(__mail_py, contE_t)
    for key in limits:
        if len(limits[key]) == 0:
            limits[key] = [-1]
    contArea = mdf.ContactAreaCONT_ELEM1D(__mail_py, contE_t)
    maxInContactLength = mdf.MaxLengthCONT_ELEM1D(__mail_py, contE_t)
    minInContactLength = mdf.MinLengthCONT_ELEM1D(__mail_py, contE_t)
    activeNodes = mdf.ActiveNodesCONT_ELEM1D(__mail_py, contE_t) # [ID, x, y]
    activeNodes_jeu = mdf.ActiveNodesFrom_JEU_ELEM1D(__mail_py, jeuE_t, maxJeu = MAX_JEU)
    activeElements = mdf.ActiveElementsCONT_ELEM1D(__mail_py, contE_t)
    xx = activeNodes[:, 1]
    yy = activeNodes[:, 2]
    ii = activeNodes[:, 0]
    xx_jeu = activeNodes_jeu[:, 1]
    yy_jeu = activeNodes_jeu[:, 2]
    ii_jeu = activeNodes_jeu[:, 0]
    table = CREA_TABLE(LISTE = (_F(PARA = 'CONT_AREA',
                                   LISTE_R = (contArea)),
                                _F(PARA = 'XLIM',
                                   LISTE_R = [xx.min(), xx.max()]),
                                _F(PARA = 'YLIM',
                                   LISTE_R = [yy.min(), yy.max()]),
                                _F(PARA = 'MIN_ELEM_LENGTH',
                                   LISTE_R = (minInContactLength)),
                                _F(PARA = 'MAX_ELEM_LENGTH',
                                   LISTE_R = (maxInContactLength)),
                                _F(PARA = 'ACTIVE_NODES',
                                    LISTE_R = (ii)),
                                _F(PARA = 'ACTIVE_NODES_JEU',
                                    LISTE_R = (ii_jeu)),
                                _F(PARA = 'ACTIVE_ELEMENTS',
                                    LISTE_I = (activeElements[:])),
                                _F(PARA = 'SHARED_NODES',
                                    LISTE_I = (limits["shared_nodes"])),
                                _F(PARA = 'SINGLE_NODES',
                                    LISTE_I = (limits["single_nodes"])),
                                _F(PARA = 'IN_ELEMENTS',
                                    LISTE_I = (limits["in_elements"])),
                                _F(PARA = 'OUT_ELEMENTS',
                                    LISTE_I = (limits["out_elements"])),
                                ))
    # }}}
    return ier
# Macro heading {{{
CONTACT_AREA_CONT_ELEM = MACRO(
        nom = 'CONTACT_AREA_CONT_ELEM',
        op = contact_area_cont_elem,
        sd_prod = table_sdaster,
        docu = "Macro to compute the contact area.",
        RESU = SIMP(statut = "o", typ = evol_noli),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        INST = SIMP(statut = "d", typ = "R", defaut = 1.0),
        MAX_JEU = SIMP(statut = "d", typ = "R", defaut = 0.0001),
)
# }}}
# }}}

# Contact area from resultat CONT_NOEU {{{
def contact_area_cont_noeu(self, RESU, MESH, INST, GROUP_MA):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("table", self.sd)
    # Add required code_aster functions
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    CREA_TABLE = self.get_cmd("CREA_TABLE")
    _F = self.get_cmd("_F")
    # }}}
    # Macro body {{{
    # Make mesh 2d
    __mail_py = MAIL_PY()
    __mail_py.FromAster(MESH)
    # Extract contact results
    __contN = CREA_CHAMP(TYPE_CHAM = 'NOEU_INFC_R',
                         OPERATION = 'EXTR',
                         RESULTAT = RESU,
                         INST = INST,
                         NOM_CHAM = 'CONT_NOEU')
    contE_jeu = __contN.EXTR_COMP('JEU', [], topo = 1)
    activeNodes = mdf.ActiveNodesJEU(contE_jeu)
    activeElements = mdf.ActiveElementsFromActiveNodes(__mail_py,
            activeNodes, GROUP_MA)
    contArea = mdf.LengthOfArrayOf1DElements(__mail_py, activeElements)
    maxInContactLength = mdf.MaxLengthOfArrayOf1DElements(__mail_py,
            activeElements)
    table = CREA_TABLE(LISTE = (_F(PARA = 'CONT_AREA',
                                   LISTE_R = (contArea)),
                                _F(PARA = 'MAX_ELEM_LENGTH',
                                   LISTE_R = (maxInContactLength)),
                                _F(PARA = 'ACTIVE_NODES',
                                   LISTE_R = (activeNodes))))
    # }}}
    return ier
# Macro heading {{{
CONTACT_AREA_CONT_NOEU = MACRO(
        nom = 'CONTACT_AREA_CONT_NOEU',
        op = contact_area_cont_noeu,
        sd_prod = table_sdaster,
        docu = "Macro to compute the contact area.",
        RESU = SIMP(statut = "o", typ = evol_noli),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        INST = SIMP(statut = "d", typ = "R", defaut = 1.0),
        GROUP_MA = SIMP(statut='o',typ=grma,validators=NoRepeat(),max='**'),
)
# }}}
# }}}

# Save contour {{{
def save_contour(self, UNITE, FILE_NAME, GROUP_MA, INST, RESU):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Add required code_aster functions
    DEFI_FICHIER = self.get_cmd("DEFI_FICHIER")
    IMPR_RESU = self.get_cmd("IMPR_RESU")
    _F = self.get_cmd("_F")
    # }}}
    # Macro body {{{
    DEFI_FICHIER(ACTION = 'ASSOCIER',
                 FICHIER = FILE_NAME,
                 UNITE = UNITE)
    if INST < 0.0:
        IMPR_RESU(RESU=_F(RESULTAT = RESU,
                          NOM_CHAM = 'DEPL',
                          GROUP_MA = GROUP_MA),
                  FORMAT = 'GMSH',
                  UNITE = UNITE)
    else:
        IMPR_RESU(RESU=_F(RESULTAT = RESU,
                          NOM_CHAM = 'DEPL',
                          INST = INST,
                          GROUP_MA = GROUP_MA),
                  FORMAT = 'GMSH',
                  UNITE = UNITE)
    DEFI_FICHIER(ACTION = 'LIBERER', UNITE = UNITE)
    # }}}
    return ier
SAVE_CONTOUR = MACRO(
        nom = 'SAVE_CONTOUR',
        op = save_contour,
        docu = "Macro to save a contour in gmsh result format.",
        UNITE=SIMP(statut='d', typ='I', defaut=101),
        FILE_NAME = SIMP(statut = "o", typ = 'TXM'),
        GROUP_MA = SIMP(statut="o", typ='TXM'),
        RESU = SIMP(statut = "o", typ = evol_noli),
        INST = SIMP(statut = "d", typ = "R", defaut = -1.0),
)
# }}}

# Mesh with connected surfaces {{{
def mesh_with_connected_surfaces(self, MESH, GROUP_NO_MAIT,
        GROUP_NO_ESCL, dummFILE, dummUNITE, GROUP_NAME):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("mesh", self.sd)
    # Add required code_aster functions
    LIRE_MAILLAGE = self.get_cmd("LIRE_MAILLAGE")
    DEFI_FICHIER = self.get_cmd("DEFI_FICHIER")
    IMPR_RESU = self.get_cmd("IMPR_RESU")
    _F = self.get_cmd("_F")
    # }}}
    # Macro body {{{
    # Test GROUP_NAME
    if len(GROUP_NAME) > 8:
        raise("Error: the length of GROUP_NAME cannot be greater than 8.")
    # Save mesh in dummFILE
    DEFI_FICHIER(ACTION = 'ASSOCIER',
                 FICHIER = dummFILE,
                 UNITE = dummUNITE)
    IMPR_RESU(UNITE = dummUNITE,
              RESU = _F(MAILLAGE = MESH),
              FORMAT = 'ASTER')
    DEFI_FICHIER(ACTION = 'LIBERER', UNITE = dummUNITE)
    # Read dummFile and create AsterMesh object
    mail = mm.AsterMesh(dummFILE)
    # Create the new elements
    mail.Add_elements_to_connect_surfaces(GROUP_NO_MAIT, GROUP_NO_ESCL,
            GROUP_NAME)
    # Write and load the new mesh
    mail.Write(dummFILE)
    DEFI_FICHIER(ACTION = 'ASSOCIER',
                 FICHIER = dummFILE,
                 UNITE = dummUNITE)
    mesh = LIRE_MAILLAGE(FORMAT = 'ASTER',
                         UNITE = dummUNITE)
    DEFI_FICHIER(ACTION = 'LIBERER', UNITE = dummUNITE)
    # }}}
    return ier
# Macro heading {{{
MESH_WITH_CONNECTED_SURFACES = MACRO(
        nom = "MESH_WITH_CONNECTED_SURFACES",
        op = mesh_with_connected_surfaces,
        sd_prod = maillage_sdaster,
        docu = "Macro command that returns a mesh with elements connecting two surfaces",
        reentrant = "n",
        fr = "",
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        GROUP_NO_MAIT = SIMP(statut="o", typ='TXM'),
        GROUP_NO_ESCL = SIMP(statut="o", typ='TXM'),
        GROUP_NAME = SIMP(statut="o", typ='TXM'),
        dummFILE = SIMP(statut="d", typ='TXM', defaut = 'dummMesh.mail'),
        dummUNITE = SIMP(statut = "d", typ = "I", defaut = 99),
        )
# }}}
# }}}

# Mesh with connected sub surfaces {{{
def mesh_with_connected_subsurfaces(self, MESH, GROUP_NO_MAIT,
        GROUP_NO_ESCL, dummFILE, dummUNITE, GROUP_NAME, NAME_SUB_MAIT,
        NAME_SUB_ESCL, SUBSURFACE, ONE_ELEMENT_PER_ESCL_NODE, XLIM,
        AUXMESH):
    # Description {{{
    """
    Add SEG2 elements from GROUP_NO_ESCL to GROUP_NO_MAIT.
    SUBSURFACE
        = 1: create a node set of the neighbouring surface nodes
        = 2: use the surface nodes.
    ONE_ELEMENT_PER_ESCL_NODE
        = 1: create one element for each slave node (to the closest
             master node)
        = 2: connect all slave nodes to all master nodes.
    """
    # }}}
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("mesh", self.sd)
    # Add required code_aster functions
    LIRE_MAILLAGE = self.get_cmd("LIRE_MAILLAGE")
    DEFI_FICHIER = self.get_cmd("DEFI_FICHIER")
    IMPR_RESU = self.get_cmd("IMPR_RESU")
    _F = self.get_cmd("_F")
    # }}}
    # Macro body {{{
    # Test GROUP_NAME
    if len(GROUP_NAME) > 8:
        raise("Error: the length of GROUP_NAME cannot be greater than 8.")
    # Save mesh in dummFILE
    if AUXMESH:
        auxFILE = 'auxMesh.mail'
        DEFI_FICHIER(ACTION = 'ASSOCIER',
                     FICHIER = auxFILE,
                     UNITE = dummUNITE)
        IMPR_RESU(UNITE = dummUNITE,
                  RESU = _F(MAILLAGE = AUXMESH),
                  FORMAT = 'ASTER')
        DEFI_FICHIER(ACTION = 'LIBERER', UNITE = dummUNITE)
        auxmail = mm.AsterMesh(auxFILE)
    DEFI_FICHIER(ACTION = 'ASSOCIER',
                 FICHIER = dummFILE,
                 UNITE = dummUNITE)
    IMPR_RESU(UNITE = dummUNITE,
              RESU = _F(MAILLAGE = MESH),
              FORMAT = 'ASTER')
    DEFI_FICHIER(ACTION = 'LIBERER', UNITE = dummUNITE)
    # Read dummFile and create AsterMesh object
    mail = mm.AsterMesh(dummFILE)
    if SUBSURFACE == 1:
        # Create sub surfaces
        mail.Add_neighbour_node_set(GROUP_NO_MAIT, NAME_SUB_MAIT)
        mail.Add_neighbour_node_set(GROUP_NO_ESCL, NAME_SUB_ESCL)
        # Create the new elements
        mail.Add_elements_to_connect_surfaces(NAME_SUB_MAIT, NAME_SUB_ESCL,
                GROUP_NAME, xlim = XLIM,
                oneElementPerSlaveNode = ONE_ELEMENT_PER_ESCL_NODE)
    else:
        # Create the new elements
        if AUXMESH:
            movedMesh = mm.AsterMesh(dummFILE)
            for nodeId in movedMesh.nodes:
                movedMesh._nodes[nodeId] = auxmail.nodes[nodeId]
            movedMesh.Add_elements_to_connect_surfaces(GROUP_NO_MAIT,
                    GROUP_NO_ESCL, GROUP_NAME, xlim = XLIM,
                    oneElementPerSlaveNode = ONE_ELEMENT_PER_ESCL_NODE,
                    rewrite = True)
            auxeles = movedMesh.elementSets[GROUP_NAME]
            eleList = []
            for eleId in auxeles:
                nods = movedMesh.elements['SEG2'][eleId]
                eleList.append(mail.add_element('SEG2', nods))
            mail._elementSets[GROUP_NAME] = eleList
        else:
            mail.Add_elements_to_connect_surfaces(GROUP_NO_MAIT, GROUP_NO_ESCL,
                    GROUP_NAME, xlim = XLIM,
                    oneElementPerSlaveNode = ONE_ELEMENT_PER_ESCL_NODE)
    # Write and load the new mesh
    mail.Write(dummFILE)
    DEFI_FICHIER(ACTION = 'ASSOCIER',
                 FICHIER = dummFILE,
                 UNITE = dummUNITE)
    mesh = LIRE_MAILLAGE(FORMAT = 'ASTER',
                         UNITE = dummUNITE)
    DEFI_FICHIER(ACTION = 'LIBERER', UNITE = dummUNITE)
    # }}}
    return ier
# Macro heading {{{
MESH_WITH_CONNECTED_SUBSURFACES = MACRO(
        nom = "MESH_WITH_CONNECTED_SUBSURFACES",
        op = mesh_with_connected_subsurfaces,
        sd_prod = maillage_sdaster,
        docu = "Macro command that returns a mesh with elements connecting two surfaces",
        reentrant = "n",
        fr = "",
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        AUXMESH = SIMP(statut = "f", typ = maillage_sdaster),
        GROUP_NO_MAIT = SIMP(statut="o", typ='TXM'),
        GROUP_NO_ESCL = SIMP(statut="o", typ='TXM'),
        GROUP_NAME = SIMP(statut="o", typ='TXM'),
        NAME_SUB_MAIT = SIMP(statut="o", typ='TXM'),
        NAME_SUB_ESCL = SIMP(statut="o", typ='TXM'),
        dummFILE = SIMP(statut="d", typ='TXM', defaut = 'dummMesh.mail'),
        dummUNITE = SIMP(statut = "d", typ = "I", defaut = 99),
        SUBSURFACE = SIMP(statut='d', typ='I', defaut= 1, into=( 1 , 2)),
        XLIM = SIMP(statut='d', typ='R', defaut= [], max ='**'),
        ONE_ELEMENT_PER_ESCL_NODE = SIMP(statut='d', typ='I', defaut= 1,
            into=( 1 , 2)),
        )
# }}}
# }}}

# Compute contact pressure from formula {{{
def calc_pression_from_formula(self, RESU, NON_DEFORMED_MESH, MODELE,
        GROUP_MA, INST, DIMS):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("pression", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    CREA_TABLE = self.get_cmd("CREA_TABLE")
    MODI_MAILLAGE = self.get_cmd("MODI_MAILLAGE")
    FORMULE = self.get_cmd("FORMULE")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    # }}}
    # Macro body {{{
    # Get stress tensor __stress {{{
    # Compute the stress field as a result
    CALC_CHAMP(reuse = RESU,
               CONTRAINTE = ('SIEF_NOEU'),
               INST = INST,
               RESULTAT = RESU)
    # Extract the result and make a field
    __stress = CREA_CHAMP(TYPE_CHAM = 'NOEU_SIEF_R',
                          OPERATION = 'EXTR',
                          RESULTAT = RESU,
                          NOM_CHAM = 'SIEF_NOEU',
                          INST = INST)
    # }}}
    # Compute normal field __normal {{{
    # Compute the displacement field that generates the deformed state
    __depl = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                        OPERATION = 'EXTR',
                        RESULTAT = RESU,
                        NOM_CHAM = 'DEPL',
                        INST = INST)
    # Compute -__depl to go from the deformed state to the undeformed one
    __mdepl = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                         OPERATION = 'COMB',
                         COMB = _F(CHAM_GD = __depl,
                                COEF_R = -1.0))
    # Deform the mesh
    MODI_MAILLAGE(reuse = NON_DEFORMED_MESH,
                  MAILLAGE = NON_DEFORMED_MESH,
                  DEFORME = _F(OPTION = 'TRAN',
                               DEPL = __depl,),);
    # Extract the normal field at the deformed state
    __normal = CREA_CHAMP(TYPE_CHAM = 'NOEU_GEOM_R',
                          OPERATION = 'NORMALE',
                          MODELE = MODELE,
                          GROUP_MA = GROUP_MA)
    # Reverse the deformation: get back to the undeformed state
    MODI_MAILLAGE(reuse = NON_DEFORMED_MESH,
                  MAILLAGE = NON_DEFORMED_MESH,
                  DEFORME = _F(OPTION = 'TRAN',
                               DEPL = __mdepl,),);
    # }}}
    # Set up contact-pressure function field __fun_pc {{{
    # Formulas for normal and tangent pressure: X and Y are the x- and
    # y-component of the normal vector
    # Functions to compute the surface traction components by Einstein notation {{{
    # Normal traction {{{
    def NormalTraction(SIXX, SIYY, SIXY, X, Y,
            SIZZ = 0.0, SIXZ = 0.0, SIYZ = 0.0, Z = 0.0):
        # Make tensors depending on each case
        if DIMS == 2:
            # Make tensors
            sig_tensor = np.array([[SIXX, SIXY], [SIXY, SIYY]])
            n_vector = np.array([X, Y])
        elif DIMS == 3:
            # Make tensors
            sig_tensor = np.array([[SIXX, SIXY, SIXZ],
                                   [SIXY, SIYY, SIYZ],
                                   [SIXZ, SIYZ, SIZZ]])
            n_vector = np.array([X, Y, Z])
        # Compute (sigma \dot n) \dot n
        val = 0.0
        for k1 in range(DIMS):
            for k2 in range(DIMS):
                val += sig_tensor[k1, k2]*n_vector[k1]*n_vector[k2]
        return -val
    # }}}
    # Tangential traction {{{
    def TangentialTraction(SIXX, SIYY, SIXY, X, Y,
            SIZZ = 0.0, SIXZ = 0.0, SIYZ = 0.0, Z = 0.0):
        # Make tensors depending on each case
        if DIMS == 2:
            # Make tensors
            sig_tensor = np.array([[SIXX, SIXY], [SIXY, SIYY]])
            n_vector = np.array([X, Y])
        elif DIMS == 3:
            # Make tensors
            sig_tensor = np.array([[SIXX, SIXY, SIXZ],
                                   [SIXY, SIYY, SIYZ],
                                   [SIXZ, SIYZ, SIZZ]])
            n_vector = np.array([X, Y, Z])
        # Compute (sigma \dot n) \dot (sigma \dot n)
        val = 0.0
        for k1 in range(DIMS):
            for k2 in range(DIMS):
                for k3 in range(DIMS):
                    val += sig_tensor[k1, k2]*n_vector[k2]*sig_tensor[k1, k3]*n_vector[k3]
        # Compute square normal traction
        pval = NormalTraction(SIXX, SIYY, SIXY, X, Y, SIZZ, SIXZ, SIYZ,
                Z)
        return abs(val - pval**2.0)**0.5
    # }}}
    # }}}
    if DIMS == 2:
        value = "func(SIXX = SIXX, SIYY = SIYY, SIXY = SIXY, X = X, Y = Y)"
        f_params = ['SIXX', 'SIYY', 'SIXY', 'X', 'Y']
    elif DIMS == 3:
        value = "func(SIXX = SIXX, SIYY = SIYY, SIXY = SIXY, X = X, Y = Y, SIZZ = SIZZ, SIXZ = SIXZ, SIYZ = SIYZ, Z = Z)"
        f_params = ['SIXX', 'SIYY', 'SIZZ', 'SIXY', 'SIXZ', 'SIYZ', 'X', 'Y', 'Z']
    __f_norm = FORMULE(VALE = value,
                       NOM_PARA = f_params,
                       func = NormalTraction)
    __f_tang = FORMULE(VALE = value,
                       NOM_PARA = f_params,
                       func = TangentialTraction)
    # Set up the contact pressure function field
    __fun_pc = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_F',
                          OPERATION = 'AFFE',
                          MAILLAGE = NON_DEFORMED_MESH,
                          AFFE=(_F(GROUP_MA = GROUP_MA,
                                  NOM_CMP = 'X1',
                                  VALE_F = __f_norm,),
                                _F(GROUP_MA = GROUP_MA,
                                  NOM_CMP = 'X2',
                                  VALE_F = __f_tang,)))
    # }}}
    # Compute the contact pressure
    pression = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                          OPERATION = 'EVAL',
                          CHAM_F = __fun_pc,
                          CHAM_PARA = (__normal, __stress))
    # }}}
    return ier
# Macro heading {{{
CALC_PRESSION_FROM_FORMULA = MACRO(
        nom = 'CALC_PRESSION_FROM_FORMULA',
        op = calc_pression_from_formula,
        sd_prod = cham_no_sdaster,
        RESU = SIMP(statut = "o", typ = evol_noli),
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        INST = SIMP(statut = "d", typ = "R", defaut = 1.0),
        DIMS = SIMP(statut = "d", typ = "I", defaut = 2, into = (2, 3)),
        GROUP_MA = SIMP(statut='o',typ=grma,validators=NoRepeat(),max='**'),
        NON_DEFORMED_MESH = SIMP(statut = "o", typ = maillage_sdaster),
        )
# }}}
# }}}

# Compute hydrostatic and shear stress {{{
def compute_hyd_shr_stress(self, RESU, TYPE, INST):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("stress", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    FORMULE = self.get_cmd("FORMULE")
    # }}}
    # Macro body {{{
    # Set formulas {{{
    if TYPE == 'D_PLAN':
        __Shyd_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ'),
                         VALE = 'Shyd(SIXX, SIYY, SIZZ)',
                         Shyd = mdf.Shyd_pln_strain)
        __Svon_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ', 'SIXY'),
                         VALE = 'Svon(SIXX, SIYY, SIZZ, SIXY)',
                         Svon = mdf.Svmis_pln_strain)
    elif TYPE == 'C_PLAN':
        __Shyd_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY'),
                         VALE = 'Shyd(SIXX, SIYY)',
                         Shyd = mdf.Shyd_pln_stress)
        __Svon_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIXY'),
                         VALE = 'Svon(SIXX, SIYY, SIXY)',
                         Svon = mdf.Svmis_pln_stress)
    else:
        raise('Error: unavailable TYPE in compute_hyd_shr_stress.')
    # }}}
    # Compute formulas in result format
    __resu = CALC_CHAMP(CONTRAINTE = ('SIGM_NOEU', 'SIGM_ELGA', 'SIEF_NOEU'),
                        CRITERES = ('SIEQ_NOEU', 'SIEQ_ELGA'),
                        INST = INST,
                        RESULTAT = RESU)
    CALC_CHAMP(CHAM_UTIL = _F(FORMULE = (__Shyd_f, __Svon_f),
                              NOM_CHAM = 'SIGM_NOEU',
                              NUME_CHAM_RESU = 1),
               INST = INST,
               reuse = __resu,
               RESULTAT = __resu)
    # Extract results in field format
    stress = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                        OPERATION = 'EXTR',
                        RESULTAT = __resu,
                        NOM_CHAM = 'UT01_NOEU')
    # }}}
    return ier
# Macro body {{{
COMPUTE_HYD_SHR_STRESS = MACRO(
        nom = 'COMPUTE_HYD_SHR_STRESS',
        op = compute_hyd_shr_stress,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute hydrostatic and shear stress fields",
        RESU = SIMP(statut = "o", typ = resultat_sdaster),
        INST = SIMP(statut = "f", typ = "R"),
        TYPE = SIMP(statut = "o", typ = 'TXM'),
        )
# }}}
# }}}

# Compute strain tensor and deviatoric strain tensor invariants {{{
def compute_strain_invariants(self, RESU, TYPE, INST):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("invs", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    FORMULE = self.get_cmd("FORMULE")
    # }}}
    # Macro body {{{
    # Set formulas {{{
    def firstInv(m11, m22, m33, m12, m13 = 0.0, m23 = 0.0):
        mat = np.array([[m11, m12, m13],
                        [m12, m22, m23],
                        [m13, m23, m33]])
        return np.trace(mat)
    def secondInv(m11, m22, m33, m12, m13 = 0.0, m23 = 0.0):
        mat = np.array([[m11, m12, m13],
                        [m12, m22, m23],
                        [m13, m23, m33]])
        mat2 = np.dot(mat, mat)
        tr2Mat = np.trace(mat)**2.0
        trMat2 = np.trace(mat2)
        return 0.5*(tr2Mat - trMat2)
    def thirdInv(m11, m22, m33, m12, m13 = 0.0, m23 = 0.0):
        mat = np.array([[m11, m12, m13],
                        [m12, m22, m23],
                        [m13, m23, m33]])
        return np.linalg.det(mat)
    def devSecondInv(m11, m22, m33, m12, m13 = 0.0, m23 = 0.0):
        mat = np.array([[m11, m12, m13],
                        [m12, m22, m23],
                        [m13, m23, m33]])
        hyd = np.trace(mat)/3.0
        mat = mat - hyd*np.identity(3)
        mat2 = np.dot(mat, mat)
        trMat2 = np.trace(mat2)
        # Return the negative second principal invariant
        return 0.5*trMat2
    def devThirdInv(m11, m22, m33, m12, m13 = 0.0, m23 = 0.0):
        mat = np.array([[m11, m12, m13],
                        [m12, m22, m23],
                        [m13, m23, m33]])
        hyd = np.trace(mat)/3.0
        mat = mat - hyd*np.identity(3)
        return np.linalg.det(mat)
    if TYPE == 'D_PLAN':
        __nom_para = ('EPXX', 'EPYY', 'EPZZ', 'EPXY')
        __input = "F(EPXX, EPYY, EPZZ, EPXY)"
    elif TYPE == 'C_PLAN':
        __nom_para = ('EPXX', 'EPYY', 'EPXY')
        __input = "F(EPXX, EPYY, 0.0, EPXY)"
    else:
        raise('Error: unavailable TYPE in compute_hyd_shr_stress.')
    __i_1_fo = FORMULE(NOM_PARA = __nom_para,
                       VALE = __input,
                       F = firstInv)
    __i_2_fo = FORMULE(NOM_PARA = __nom_para,
                       VALE = __input,
                       F = secondInv)
    __i_3_fo = FORMULE(NOM_PARA = __nom_para,
                       VALE = __input,
                       F = thirdInv)
    __d_2_fo = FORMULE(NOM_PARA = __nom_para,
                       VALE = __input,
                       F = devSecondInv)
    __d_3_fo = FORMULE(NOM_PARA = __nom_para,
                       VALE = __input,
                       F = devThirdInv)
    # }}}
    # Compute formulas in result format
    __resu = CALC_CHAMP(DEFORMATION = ('EPSI_NOEU'),
                        INST = INST,
                        RESULTAT = RESU)
    CALC_CHAMP(CHAM_UTIL = _F(FORMULE = (__i_1_fo,
                                         __i_2_fo,
                                         __i_3_fo,
                                         __d_2_fo,
                                         __d_3_fo),
                              NOM_CHAM = 'EPSI_NOEU',
                              NUME_CHAM_RESU = 1),
               INST = INST,
               reuse = __resu,
               RESULTAT = __resu)
    # Extract results in field format
    invs = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                        OPERATION = 'EXTR',
                        RESULTAT = __resu,
                        NOM_CHAM = 'UT01_NOEU')
    # }}}
    return ier
# Macro body {{{
COMPUTE_STRAIN_INVARIANTS = MACRO(
        nom = 'COMPUTE_STRAIN_INVARIANTS',
        op = compute_strain_invariants,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute the strain invariants",
        RESU = SIMP(statut = "o", typ = resultat_sdaster),
        INST = SIMP(statut = "f", typ = "R"),
        TYPE = SIMP(statut = "o", typ = 'TXM'),
        )
# }}}
# }}}

# Morphogenesis growth stress with function {{{
def morphogenesis_growth_with_function(self, TYPE, RESU, INST, GROUP_MA,
        MODELE, GROWTH_FUNCTION):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("growth_t", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    FORMULE = self.get_cmd("FORMULE")
    # }}}
    # Macro body {{{
    # Compute necessary stress fields {{{
    CALC_CHAMP(reuse = RESU,
               CONTRAINTE = ('SIGM_NOEU', 'SIGM_ELGA',
                             'SIEF_NOEU'),
               CRITERES = ('SIEQ_ELGA'),
               INST = INST,
               RESULTAT = RESU)
    # }}}
    # Compute hydrostatic and shear stress {{{
    # Formulas
    if TYPE == 'D_PLAN':
        __Shyd_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ'),
                         VALE = 'Shyd(SIXX, SIYY, SIZZ)',
                         Shyd = mdf.Shyd_pln_strain)
        __Svon_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ', 'SIXY'),
                         VALE = 'Svon(SIXX, SIYY, SIZZ, SIXY)',
                         Svon = mdf.Svmis_pln_strain)
    elif TYPE == 'C_PLAN':
        __Shyd_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY'),
                         VALE = 'Shyd(SIXX, SIYY)',
                         Shyd = mdf.Shyd_pln_stress)
        __Svon_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIXY'),
                         VALE = 'Svon(SIXX, SIYY, SIXY)',
                         Svon = mdf.Svmis_pln_stress)
    else:
        raise('Error: unavailable TYPE in morphogenesis_growth_with_function.')
    # Stresses as result
    __resu = CALC_CHAMP(CHAM_UTIL = _F(FORMULE = (__Shyd_f, __Svon_f),
                                       NOM_CHAM = 'SIGM_NOEU',
                                       NUME_CHAM_RESU = 1),
                        RESULTAT = RESU,
                        INST = INST)
    # Stresses as fields
    __stress = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                          OPERATION = 'EXTR',
                          RESULTAT = __resu,
                          INST = INST,
                          NOM_CHAM = 'UT01_NOEU')
    # }}}
    # Compute growth stress tensor {{{
    # Function field
    __grow_f = CREA_CHAMP(OPERATION = 'AFFE',
                          TYPE_CHAM = 'NOEU_NEUT_F',
                          MODELE = MODELE,
                          AFFE = _F(GROUP_MA = GROUP_MA,
                                    NOM_CMP = 'X1',
                                    VALE_F = GROWTH_FUNCTION))
    # Evaluate function
    __grow = CREA_CHAMP(OPERATION = 'EVAL',
                        TYPE_CHAM = 'NOEU_NEUT_R',
                        CHAM_F = __grow_f,
                        CHAM_PARA = __stress)
    # Create tensor field
    FFF = [_F(GROUP_MA = GROUP_MA,
             CHAM_GD = __grow,
             NOM_CMP = ('X1'),
             NOM_CMP_RESU = ('SIXX')),
           _F(GROUP_MA = GROUP_MA,
             CHAM_GD = __grow,
             NOM_CMP = ('X1'),
             NOM_CMP_RESU = ('SIYY'))]
    if TYPE == 'D_PLAN':
        FFF.append(_F(GROUP_MA = GROUP_MA,
                      CHAM_GD = __grow,
                      NOM_CMP = ('X1'),
                      NOM_CMP_RESU = ('SIZZ')))
    elif TYPE == 'C_PLAN':
        pass
    else:
        raise('Error: unavailable TYPE in morphogenesis_growth_with_function.')
    growth_t = CREA_CHAMP(OPERATION = 'ASSE',
                          TYPE_CHAM = 'NOEU_SIEF_R',
                          MODELE = MODELE,
                          ASSE = FFF)
    # }}}
    # }}}
    return ier
# Macro heading {{{
MORPHOGENESIS_GROWTH_WITH_FUNCTION = MACRO(
        nom = 'MORPHOGENESIS_GROWTH_WITH_FUNCTION',
        op = morphogenesis_growth_with_function,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute the morphogenesis-inspired growth stress tensor.",
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        GROUP_MA = SIMP(statut = "o", typ = 'TXM'),
        INST = SIMP(statut = "f", typ = "R"),
        RESU = SIMP(statut = "o", typ = resultat_sdaster),
        TYPE = SIMP(statut = "o", typ = 'TXM'),
        GROWTH_FUNCTION = SIMP(statut='o', typ=formule, max='**'),
        )
# }}}
# }}}

# Compute mesh morphogenesis displacement with function {{{
def compute_mesh_morphogenesis_displacement_with_function(self, RESU,
        INST, GROUP_MA, MODELE, EXCIT, CHAM_MATER, REF_FIXED_DX,
        REF_FIXED_DY, CARA_ELEM, GROWTH_FUNCTION, TYPE):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("disp", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    AFFE_CHAR_MECA = self.get_cmd("AFFE_CHAR_MECA")
    MECA_STATIQUE = self.get_cmd("MECA_STATIQUE")
    FORMULE = self.get_cmd("FORMULE")
    # }}}
    # Macro body {{{
    # Tests group inputs {{{
    if isinstance(REF_FIXED_DX, str):
        fixedDX = REF_FIXED_DX
    else:
        if REF_FIXED_DX == None:
            fixedDX = None
        else:
            if len(REF_FIXED_DX) == 1:
                fixedDX = REF_FIXED_DX[0]
            elif len(REF_FIXED_DX) == 0:
                fixedDX = None
            else:
                raise("Error: REF_FIXED_DX includes several groups.")
    if isinstance(REF_FIXED_DY, str):
        fixedDY = REF_FIXED_DY
    else:
        if REF_FIXED_DY == None:
            fixedDY = None
        else:
            if len(REF_FIXED_DY) == 1:
                fixedDY = REF_FIXED_DY[0]
            elif len(REF_FIXED_DY) == 0:
                fixedDY = None
            else:
                raise("Error: REF_FIXED_DY includes several groups.")
    # }}}
    # Compute growth tensor {{{
    __grw_t = MORPHOGENESIS_GROWTH_WITH_FUNCTION(
            TYPE = TYPE,
            RESU = RESU,
            INST = INST,
            GROUP_MA = GROUP_MA,
            MODELE = MODELE,
            GROWTH_FUNCTION = GROWTH_FUNCTION
            )
    __grw_g = CREA_CHAMP(OPERATION = 'DISC',
                         TYPE_CHAM = 'ELGA_SIEF_R',
                         MODELE = MODELE,
                         CHAM_GD = __grw_t,
                         PROL_ZERO = 'OUI')
    # }}}
    # Make growth load {{{
    __preStr = AFFE_CHAR_MECA(MODELE = MODELE,
                              PRE_SIGM = _F(SIGM = __grw_g))
    loads = []
    for load in EXCIT:
        loads.append(load)
    loads.append(_F(CHARGE = __preStr))
    # }}}
    # Solve growth linear problem {{{
    __resu = MECA_STATIQUE(
            CHAM_MATER = CHAM_MATER,
            CARA_ELEM = CARA_ELEM,
            EXCIT = loads,
            MODELE = MODELE
            )
    # }}}
    # Compute growth displacement {{{
    # Formula
    __disp1 = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                         OPERATION = 'EXTR',
                         RESULTAT = __resu,
                         NOM_CHAM = 'DEPL')
    # Get the displacement assign to the reference node, the node
    # where the displacement should be zero
    # In the x-direction
    if fixedDX == None:
        refDx = 0.0
    else:
        refDx = __disp1.EXTR_COMP('DX', [fixedDX]).valeurs
        if len(refDx) == 1:
            refDx = refDx[0]
        else:
            raise("Error: the reference group DX includes several nodes.")
    # In the y-direction
    if fixedDY == None:
        refDy = 0.0
    else:
        refDy = __disp1.EXTR_COMP('DY', [fixedDY]).valeurs
        if len(refDy) == 1:
            refDy = refDy[0]
        else:
            raise("Error: the reference group DY includes several nodes.")
    # Set formulas
    __dispxf = FORMULE(NOM_PARA = ('DX'), VALE = 'DX - refDx', refDx = refDx)
    __dispyf = FORMULE(NOM_PARA = ('DY'), VALE = 'DY - refDy', refDy = refDy)
    __dispze = FORMULE(NOM_PARA = ('DX','DY'), VALE = '0.0')
    # Set function field
    __disp_f = CREA_CHAMP(OPERATION = 'AFFE',
                          TYPE_CHAM = 'NOEU_NEUT_F',
                          MODELE = MODELE,
                          AFFE = (_F(TOUT = 'OUI',
                                     NOM_CMP = ('X1', 'X2'),
                                     VALE_F = (__dispze, __dispze)),
                                  _F(GROUP_MA = GROUP_MA,
                                     NOM_CMP = ('X1', 'X2'),
                                     VALE_F = (__dispxf, __dispyf))))
    # Evaluate displacements only in GROUP_MA
    __disp2 = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                         OPERATION = 'EVAL',
                         CHAM_F = __disp_f,
                         CHAM_PARA = __disp1)
    # Convert the neutral field into a displacement one
    disp = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                      OPERATION = 'ASSE',
                      MODELE = MODELE,
                      ASSE = _F(CHAM_GD = __disp2,
                                TOUT = 'OUI',
                                NOM_CMP = ('X1', 'X2'),
                                NOM_CMP_RESU = ('DX', 'DY')))
    # }}}
    # }}}
    return ier
# Macro heading {{{
COMPUTE_MESH_MORPHOGENESIS_DISPLACEMENT_WITH_FUNCTION = MACRO(
        nom = 'COMPUTE_MESH_MORPHOGENESIS_DISPLACEMENT_WITH_FUNCTION',
        op = compute_mesh_morphogenesis_displacement_with_function,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute the morphogenesis-inspired displacement.",
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        GROUP_MA = SIMP(statut = "o", typ = 'TXM'),
        REF_FIXED_DX = SIMP(statut='f',typ=grma,validators=NoRepeat(),max='**', min = 0),
        REF_FIXED_DY = SIMP(statut='f',typ=grma,validators=NoRepeat(),max='**', min = 0),
        INST = SIMP(statut = "f", typ = "R"),
        RESU = SIMP(statut = "o", typ = resultat_sdaster),
        CHAM_MATER = SIMP(statut = 'o', typ=cham_mater),
        EXCIT = FACT(statut = 'o', max='**',
                     CHARGE = SIMP(statut = 'o', typ = (char_meca,
                         char_cine_meca)),
                     FONC_MULT = SIMP(statut = 'f', typ = (fonction_sdaster,
                         nappe_sdaster, formule)),
                     TYPE_CHARGE = SIMP(statut='f', typ='TXM', defaut = "FIXE_CSTE",
                         into=("FIXE_CSTE", "FIXE_PILO", "SUIV",
                             "SUIV_PILO", "DIDI"))),
        CARA_ELEM =SIMP(statut='f',typ=cara_elem),
        TYPE = SIMP(statut = "o", typ = 'TXM'),
        GROWTH_FUNCTION = SIMP(statut='o', typ=formule, max='**'),
        )
# }}}
# }}}

# Deviatoric stress {{{
def deviatoric_stress(self, RESU, MODELE, TYPE, INST):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("devia", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    FORMULE = self.get_cmd("FORMULE")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    # }}}
    # Macro body {{{
    __resu = CALC_CHAMP(RESULTAT = RESU,
                        CONTRAINTE = 'SIEF_NOEU')
    __stress = CREA_CHAMP(TYPE_CHAM = 'NOEU_SIEF_R',
                          INST = INST,
                          OPERATION = 'EXTR',
                          RESULTAT = __resu,
                          NOM_CHAM = 'SIEF_NOEU')
    if TYPE == 'D_PLAN':
        __hyd_fo = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ'),
                           VALE = '(SIXX + SIYY + SIZZ)/3.0')
    elif TYPE == 'C_PLAN':
        __hyd_fo = FORMULE(NOM_PARA = ('SIXX', 'SIYY'),
                           VALE = '(SIXX + SIYY)/3.0')
    else:
        raise('Error: unavailable TYPE in deviatoric_stress.')
    __zer_fo = FORMULE(NOM_PARA = ('SIXX', 'SIYY'),
                       VALE = '0.0')
    __hyd_fu = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_F',
                          OPERATION = 'AFFE',
                          MODELE = MODELE,
                          AFFE = (_F(TOUT = 'OUI',
                                     NOM_CMP = ('X1', 'X2'),
                                     VALE_F = (__hyd_fo, __zer_fo))))
    __hyd_va = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                          OPERATION = 'EVAL',
                          CHAM_F = __hyd_fu,
                          CHAM_PARA = __stress)
    __hyd_fi = CREA_CHAMP(OPERATION = 'ASSE',
                          TYPE_CHAM = 'NOEU_SIEF_R',
                          MODELE = MODELE,
                          ASSE = (_F(TOUT = 'OUI',
                                     CHAM_GD = __hyd_va,
                                     NOM_CMP = 'X1',
                                     NOM_CMP_RESU = 'SIXX'),
                                  _F(TOUT = 'OUI',
                                     CHAM_GD = __hyd_va,
                                     NOM_CMP = 'X1',
                                     NOM_CMP_RESU = 'SIYY'),
                                  _F(TOUT = 'OUI',
                                     CHAM_GD = __hyd_va,
                                     NOM_CMP = 'X1',
                                     NOM_CMP_RESU = 'SIZZ'),
                                  _F(TOUT = 'OUI',
                                     CHAM_GD = __hyd_va,
                                     NOM_CMP = 'X2',
                                     NOM_CMP_RESU = 'SIXY')))
    devia = CREA_CHAMP(OPERATION = 'COMB',
                       TYPE_CHAM = 'NOEU_SIEF_R',
                       COMB = (_F(CHAM_GD = __stress,
                                  COEF_R = 1.0),
                               _F(CHAM_GD = __hyd_fi,
                                  COEF_R = -1.0)))
    # }}}
    return ier
# Macro heading {{{
DEVIATORIC_STRESS = MACRO(
        nom = 'DEVIATORIC_STRESS',
        op = deviatoric_stress,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute the deviatoric stress tensor.",
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        INST = SIMP(statut = "f", typ = "R"),
        RESU = SIMP(statut = "o", typ = resultat_sdaster),
        TYPE = SIMP(statut = "o", typ = 'TXM'),
        )
# }}}
# }}}

# Compute equivalent uniform force {{{
def compute_equivalent_uniform_force(self, MODELE, GROUP_MA_CONT, RESU,
        MESH, INST, GROUP_NO_DX_0, GROUP_NO_DY_0):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("charge", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    AFFE_CHAR_MECA = self.get_cmd("AFFE_CHAR_MECA")
    # }}}
    # Macro body {{{
    # Get normal
    __normal = CREA_CHAMP(TYPE_CHAM = 'NOEU_GEOM_R',
                          OPERATION = 'NORMALE',
                          MODELE = MODELE,
                          GROUP_MA = (GROUP_MA_CONT))
    # Get contact area and active nodes
    __tb_co = CONTACT_AREA_CONT_ELEM(RESU = RESU,
                                      MESH = MESH,
                                      INST = INST)
    activeElements = __tb_co.EXTR_TABLE().values()["ACTIVE_ELEMENTS"]
    # Get contact pressure
    try:
        __tbcont = CONTACT_PRESSURE_ABSC_CURV(RESU = RESU,
                                         MESH = MESH,
                                         GROUP_MA = GROUP_MA_CONT,
                                         INST = INST)
        # Get pressure vector
        p_c = np.array(__tbcont.EXTR_TABLE().values()['LAGS_C'])
    except:
        __tbcont = CONTACT_PRESSURE_ABSC_CURV(RESU = RESU,
                                         MESH = MESH,
                                         GROUP_MA = GROUP_MA_CONT,
                                         MODELE = MODELE,
                                         INST = INST)
        # Get pressure vector
        p_c = np.array(__tbcont.EXTR_TABLE().values()['X1'])
    xx = np.array(__tbcont.EXTR_TABLE().values()['ABSC_CURV'])
    # Get nodal force _F
    __mail_py = MAIL_PY()
    __mail_py.FromAster(MESH)
    FFF1 = ElementEquivalentPressure_F(xx, p_c, __tb_co, __mail_py, area = 0.454015)
    # Set up displacement bc
    FFF2 = [_F(DX = 0.0, GROUP_NO = GROUP_NO_DX_0),
            _F(DY = 0.0, GROUP_NO = GROUP_NO_DY_0)]
    # Set up load
    charge =  AFFE_CHAR_MECA(
            PRES_REP = (FFF1),
            MODELE = MODELE)
    # }}}
    return ier
# Macro heading {{{
COMPUTE_EQUIVALENT_UNIFORM_FORCE = MACRO(
        nom = 'COMPUTE_EQUIVALENT_UNIFORM_FORCE',
        op = compute_equivalent_uniform_force,
        sd_prod = char_meca,
        docu = "Macro to compute a displacement field trigger by wear.",
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        GROUP_MA_CONT = SIMP(statut = "o", typ = 'TXM'),
        INST = SIMP(statut = "d", typ = "R", defaut = 1.0),
        RESU = SIMP(statut = "o", typ = evol_noli),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        GROUP_NO_DX_0 = SIMP(statut='f',typ=grma ,validators=NoRepeat(),max='**'),
        GROUP_NO_DY_0 = SIMP(statut='f',typ=grma ,validators=NoRepeat(),max='**'),
        )
# }}}
# }}}

# Compute displacement around the growth zone {{{
def compute_mesh_growth_around_contact(self, MODELE, GROUP_MA_CONT,
        RESU, MESH, INST, GROUP_NO_DX_0, GROUP_NO_DY_0, CHAM_MATER,
        CARA_ELEM, GROWTH_DISP, GROUP_MA, MAX_JEU, NORMAL_DIRECTION):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("disp", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    AFFE_CHAR_MECA = self.get_cmd("AFFE_CHAR_MECA")
    MECA_STATIQUE = self.get_cmd("MECA_STATIQUE")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    # }}}
    # Macro body {{{
    # Get active nodes and contact area
    __tbcoar = CONTACT_AREA_CONT_ELEM(RESU = RESU,
                                      MESH = MESH,
                                      INST = INST)
    contArea = __tbcoar.EXTR_TABLE().values()["CONT_AREA"][0]
    activeNodes = __tbcoar.EXTR_TABLE().values()["ACTIVE_NODES"]
    activeNodes = list(filter((None).__ne__, activeNodes))
    activeNodes = [int(val) + 1 for val in activeNodes] # + 1 to set first node as 1 and not 0
    activeNodes_jeu = __tbcoar.EXTR_TABLE().values()["ACTIVE_NODES_JEU"]
    activeNodes_jeu = list(filter((None).__ne__, activeNodes_jeu))
    activeNodes_jeu = [int(val) + 1 for val in activeNodes_jeu] # + 1 to set first node as 1 and not 0
    # Get displacement and normal at GROUP_MA_CONT
    if isinstance(GROUP_MA_CONT, str):
        lgma = [GROUP_MA_CONT]
    else:
        lgma = GROUP_MA_CONT
    disp_contact_dx = GROWTH_DISP.EXTR_COMP('DX', lgma, topo = 1)
    disp_contact_dy = GROWTH_DISP.EXTR_COMP('DY', lgma, topo = 1)
    # Get normal if NORMAL_DIRECTION = 2
    if NORMAL_DIRECTION == 2:
        __normal = CREA_CHAMP(TYPE_CHAM = 'NOEU_GEOM_R',
                              OPERATION = 'NORMALE',
                              MODELE = MODELE,
                              GROUP_MA = (GROUP_MA_CONT))
        nx = __normal.EXTR_COMP('X', lgma, topo = 1)
        ny = __normal.EXTR_COMP('Y', lgma, topo = 1)
    # Set displacement on active nodes
    __mail_py = MAIL_PY()
    __mail_py.FromAster(MESH)
    disp_FFF = []
    noeud = disp_contact_dx.noeud
    for k1 in range(len(noeud)):
        node_Id = noeud[k1]
        if node_Id in activeNodes_jeu: # Active nodes by jeu.
            node_name = __mail_py.correspondance_noeuds[node_Id - 1]
            dx_i = disp_contact_dx.valeurs[k1]
            dy_i = disp_contact_dy.valeurs[k1]
            if NORMAL_DIRECTION == 1:
                ux_i = dx_i
                uy_i = dy_i
            else:
                nx_i = nx.valeurs[k1]
                ny_i = ny.valeurs[k1]
                dot = (nx_i*dx_i + ny_i*dy_i)
                ux_i = dot*nx_i
                uy_i = dot*ny_i
            disp_FFF.append(_F(NOEUD = node_name,
                               DX = ux_i,
                               DY = uy_i))
    # Add zero displacement on GROUP_NO_DX_0, and GROUP_NO_DY_0
    disp_FFF.append(_F(DX = 0.0, GROUP_NO = GROUP_NO_DX_0))
    disp_FFF.append(_F(DY = 0.0, GROUP_NO = GROUP_NO_DY_0))
    # Compute displacement
    __load =  AFFE_CHAR_MECA(DDL_IMPO = disp_FFF,
                             MODELE = MODELE)
    __resu = MECA_STATIQUE(CHAM_MATER = CHAM_MATER,
                           EXCIT = (_F(CHARGE = __load)),
                           CARA_ELEM = CARA_ELEM,
                           MODELE = MODELE)
    disp = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                      OPERATION = 'EXTR',
                      RESULTAT = __resu,
                      NOM_CHAM = 'DEPL')
    # }}}
    return ier
# Macro heading {{{
COMPUTE_MESH_GROWTH_AROUND_CONTACT = MACRO(
        nom = 'COMPUTE_MESH_GROWTH_AROUND_CONTACT',
        op = compute_mesh_growth_around_contact,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute a displacement field around the contact area.",
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        GROUP_MA_CONT = SIMP(statut = "o", typ = 'TXM'),
        RESU = SIMP(statut = "o", typ = evol_noli),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        INST = SIMP(statut = "d", typ = "R", defaut = 1.0),
        GROUP_NO_DX_0 = SIMP(statut='f',typ=grma ,validators=NoRepeat(),max='**'),
        GROUP_NO_DY_0 = SIMP(statut='f',typ=grma ,validators=NoRepeat(),max='**'),
        CHAM_MATER = SIMP(statut = 'o', typ=cham_mater),
        CARA_ELEM =SIMP(statut='f',typ=cara_elem),
        GROUP_MA = SIMP(statut = "o", typ = 'TXM'),
        GROWTH_DISP = SIMP(statut = "o", typ = cham_no_sdaster),
        MAX_JEU = SIMP(statut = "d", typ = "R", defaut = 0.0001),
        NORMAL_DIRECTION = SIMP(statut = 'd', typ = 'I', defaut = 1, into = (1, 2)),
        )
# }}}
# }}}

# Apply distributed total force contour {{{
def apply_2D_total_force_contour(self, MESH, MODELE,
        TOTAL_FORCE_CONTOUR, TIME):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut('char', self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    AFFE_CHAR_MECA = self.get_cmd("AFFE_CHAR_MECA")
    AFFE_CHAR_MECA_F = self.get_cmd("AFFE_CHAR_MECA_F")
    DEFI_FONCTION = self.get_cmd("DEFI_FONCTION")
    # }}}
    # Macro body {{{
    # Load mesh
    __mail_py = MAIL_PY()
    __mail_py.FromAster(MESH)
    nodes = __mail_py.cn
    connectivities = __mail_py.co
    # Get equivalent pressure for each force {{{
    f_dirs = {"FX", "FY"}
    force_contours = []
    for f_i in TOTAL_FORCE_CONTOUR:
        # Get force and moment
        f_dir_i = list(f_dirs.intersection(f_i.keys()))
        if len(f_dir_i) != 1:
            message = "The ABS_FORCE_FACE " + str(f_i) + " must have "
            message += "and only one key 'FX', 'FY' or 'FZ'."
            raise ValueError(message)
        f_dir_i = f_dir_i[0]
        M_i = f_i.get("M", 0.0)
        f_mag_i = f_i[f_dir_i]
        # Get period and phase
        factor_function = f_i.get("FACTOR_FUNCTION", "lambda t: 1.0")
        factor_parameters = f_i.get("FACTOR_PARAMETERS", {})
        factor_function = eval(factor_function)
        # Compute factor for time t (if period is well defined)
        factor_i = factor_function(TIME, **factor_parameters)
        # Get group of elements
        groupma_i = f_i["GROUP_MA"]
        # Convert string to list of strings if necessary
        if isinstance(groupma_i, str):
            groupma_i = [groupma_i]
        # Find the total length of the contour
        len_i = 0.0
        for name_j in groupma_i:
            len_i += mdf.Length1DElementGroup(__mail_py, name_j)
        # Add information
        force_contours.append({"DIR" : f_dir_i,
                               "MAG" : f_mag_i,
                               "MOM" : M_i*factor_i,
                               "LEN" : len_i,
                               "GROUP_MA" : groupma_i
                              })
    # }}}
    # Make _F() for each force contour {{{
    force_contour_Fs = []
    _fo_fu = []
    for f_i in force_contours:
        mom_i = f_i["MOM"]
        mag_i = f_i["MAG"]
        dir_i = f_i["DIR"]
        len_i = f_i["LEN"]
        groupma_i = f_i["GROUP_MA"]
        # Get slope (m) and cut point (b)
        m = -12.0*mom_i/(len_i**3.0)
        b = mag_i/len_i - m*len_i/2.0
        # Get F(0) and F(L)
        F_0 = b
        F_L = m*len_i + b
        # Get spatial limits {{{
        # Select direction of the limits (opposite to the load)
        if dir_i == "FX":
            # Get limits in the y direction
            coord_id = 1
            nom_para = 'Y'
        else:
            # Get limits in the x direction
            coord_id = 0
            nom_para = 'X'
        print("************** Warning ****************")
        print("In MACRO ABS_FORCE_FACE")
        print("The force is always defined as a function of 'X'.")
        print("***************************************")
        coord_id = 0
        nom_para = "X"
        # Get the limits
        lmax = -9999999999.9
        lmin =  9999999999.9
        for name_j in groupma_i:
            group_eles_i = __mail_py.gma[name_j]
            for eleId in group_eles_i:
                eleNods = connectivities[eleId]
                for nodeId in eleNods:
                    coord_i = nodes[nodeId][coord_id]
                    if coord_i > lmax:
                        lmax = coord_i
                    if coord_i < lmin:
                        lmin = coord_i
        # }}}
        # Set linearly distributed force
        _fo_fu.append(DEFI_FONCTION(NOM_PARA = nom_para,
                                    VALE = (lmin, F_0, lmax, F_L)))
        dict_F = {dir_i : _fo_fu[-1],
                  "GROUP_MA" : groupma_i}
        # Add _F
        force_contour_Fs.append(_F(**dict_F))
    # }}}
    char = AFFE_CHAR_MECA_F(FORCE_CONTOUR = force_contour_Fs,
                          MODELE = MODELE)
    # }}}
    return ier
# Macro heading {{{
APPLY_2D_TOTAL_FORCE_CONTOUR = MACRO(
        nom = "APPLY_2D_TOTAL_FORCE_CONTOUR",
        op = apply_2D_total_force_contour,
        sd_prod = char_meca,
        docu = "Macro that returns the application of a total force at a contour.",
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        MODELE = SIMP(statut='o', typ = modele_sdaster),
        TOTAL_FORCE_CONTOUR = SIMP(statut = 'f', typ = dict, max = '**'),
        TIME = SIMP(statut = 'd', typ = "R", defaut = 1.0),
        )
# }}}
# }}}

# Compute displacements from growth force {{{
def compute_mesh_displacement_from_growth_force(self, GRW_TEN, MODELE,
        EXCIT, CHAM_MATER, CARA_ELEM, REF_FIXED_DX, REF_FIXED_DY,
        GROUP_MA):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("disp", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    AFFE_CHAR_MECA = self.get_cmd("AFFE_CHAR_MECA")
    MECA_STATIQUE = self.get_cmd("MECA_STATIQUE")
    FORMULE = self.get_cmd("FORMULE")
    # }}}
    # Macro body {{{
    # Tests group inputs {{{
    if isinstance(REF_FIXED_DX, str):
        fixedDX = REF_FIXED_DX
    else:
        if REF_FIXED_DX == None:
            fixedDX = None
        else:
            if len(REF_FIXED_DX) == 1:
                fixedDX = REF_FIXED_DX[0]
            elif len(REF_FIXED_DX) == 0:
                fixedDX = None
            else:
                raise("Error: REF_FIXED_DX includes several groups.")
    if isinstance(REF_FIXED_DY, str):
        fixedDY = REF_FIXED_DY
    else:
        if REF_FIXED_DY == None:
            fixedDY = None
        else:
            if len(REF_FIXED_DY) == 1:
                fixedDY = REF_FIXED_DY[0]
            elif len(REF_FIXED_DY) == 0:
                fixedDY = None
            else:
                raise("Error: REF_FIXED_DY includes several groups.")
    # }}}
    # Compute growth tensor {{{
    __grw_g = CREA_CHAMP(OPERATION = 'DISC',
                         TYPE_CHAM = 'ELGA_SIEF_R',
                         MODELE = MODELE,
                         CHAM_GD = GRW_TEN,
                         PROL_ZERO = 'OUI')
    # }}}
    # Make growth load {{{
    __preStr = AFFE_CHAR_MECA(MODELE = MODELE,
                              PRE_SIGM = _F(SIGM = __grw_g))
    loads = []
    for load in EXCIT:
        loads.append(load)
    loads.append(_F(CHARGE = __preStr))
    # }}}
    # Solve growth linear problem {{{
    __resu = MECA_STATIQUE(
            CHAM_MATER = CHAM_MATER,
            CARA_ELEM = CARA_ELEM,
            EXCIT = loads,
            MODELE = MODELE
            )
    # }}}
    # Compute growth displacement {{{
    # Formula
    __disp1 = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                         OPERATION = 'EXTR',
                         RESULTAT = __resu,
                         NOM_CHAM = 'DEPL')
    # Get the displacement assign to the reference node, the node
    # where the displacement should be zero
    # In the x-direction
    if fixedDX == None:
        refDx = 0.0
    else:
        refDx = __disp1.EXTR_COMP('DX', [fixedDX]).valeurs
        if len(refDx) == 1:
            refDx = refDx[0]
        else:
            raise("Error: the reference group DX includes several nodes.")
    # In the y-direction
    if fixedDY == None:
        refDy = 0.0
    else:
        refDy = __disp1.EXTR_COMP('DY', [fixedDY]).valeurs
        if len(refDy) == 1:
            refDy = refDy[0]
        else:
            raise("Error: the reference group DY includes several nodes.")
    # Set formulas
    __dispxf = FORMULE(NOM_PARA = ('DX'), VALE = 'DX - refDx', refDx = refDx)
    __dispyf = FORMULE(NOM_PARA = ('DY'), VALE = 'DY - refDy', refDy = refDy)
    __dispze = FORMULE(NOM_PARA = ('DX','DY'), VALE = '0.0')
    # Set function field
    __disp_f = CREA_CHAMP(OPERATION = 'AFFE',
                          TYPE_CHAM = 'NOEU_NEUT_F',
                          MODELE = MODELE,
                          AFFE = (_F(TOUT = 'OUI',
                                     NOM_CMP = ('X1', 'X2'),
                                     VALE_F = (__dispze, __dispze)),
                                  _F(GROUP_MA = GROUP_MA,
                                     NOM_CMP = ('X1', 'X2'),
                                     VALE_F = (__dispxf, __dispyf))))
    # Evaluate displacements only in GROUP_MA
    __disp2 = CREA_CHAMP(TYPE_CHAM = 'NOEU_NEUT_R',
                         OPERATION = 'EVAL',
                         CHAM_F = __disp_f,
                         CHAM_PARA = __disp1)
    # Convert the neutral field into a displacement one
    disp = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                      OPERATION = 'ASSE',
                      MODELE = MODELE,
                      ASSE = _F(CHAM_GD = __disp2,
                                TOUT = 'OUI',
                                NOM_CMP = ('X1', 'X2'),
                                NOM_CMP_RESU = ('DX', 'DY')))
    # }}}
    # }}}
    return
# Macro heading {{{
COMPUTE_MESH_DISPLACEMENT_FROM_GROWTH_FORCE = MACRO(
        nom = 'COMPUTE_MESH_DISPLACEMENT_FROM_GROWTH_FORCE',
        op = compute_mesh_displacement_from_growth_force,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute the displacement due to a growth force.",
        GRW_TEN = SIMP(statut = 'o', typ = cham_no_sdaster),
        MODELE = SIMP(statut='o', typ = modele_sdaster),
        GROUP_MA = SIMP(statut = "o", typ = 'TXM'),
        REF_FIXED_DX = SIMP(statut='f',typ=grma,validators=NoRepeat(),max='**', min = 0),
        REF_FIXED_DY = SIMP(statut='f',typ=grma,validators=NoRepeat(),max='**', min = 0),
        CHAM_MATER = SIMP(statut = 'o', typ=cham_mater),
        EXCIT = FACT(statut = 'o', max='**',
                     CHARGE = SIMP(statut = 'o', typ = (char_meca,
                         char_cine_meca)),
                     FONC_MULT = SIMP(statut = 'f', typ = (fonction_sdaster,
                         nappe_sdaster, formule)),
                     TYPE_CHARGE = SIMP(statut='f', typ='TXM', defaut = "FIXE_CSTE",
                         into=("FIXE_CSTE", "FIXE_PILO", "SUIV",
                             "SUIV_PILO", "DIDI"))),
        CARA_ELEM =SIMP(statut='f',typ=cara_elem))
# }}}
# }}}

# Morphogenesis growth stress beneath contour {{{
def morphogenesis_growth_beneath_contour(self, RESU, INST, GROUP_MA,
        GROUP_MA_CONT, MODELE, MODELISATION, MESH, GFUNC, FUNC_PARAMS,
        SMOOTH_PARAMS, GEOMETRIE):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("growth_t", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    CALC_CHAMP = self.get_cmd("CALC_CHAMP")
    FORMULE = self.get_cmd("FORMULE")
    MODI_MAILLAGE = self.get_cmd("MODI_MAILLAGE")
    # }}}
    # Macro body {{{
    # Compute necessary stress fields {{{
    CALC_CHAMP(reuse = RESU,
               CONTRAINTE = ('SIGM_NOEU', 'SIGM_ELGA',
                             'SIEF_NOEU'),
               CRITERES = ('SIEQ_ELGA'),
               INST = INST,
               RESULTAT = RESU)
    __sigma = CREA_CHAMP(TYPE_CHAM = 'NOEU_SIEF_R',
                         OPERATION = 'EXTR',
                         RESULTAT = RESU,
                         NOM_CHAM = 'SIEF_NOEU',
                         INST = INST)
    # }}}
    growth_t = MG_GROWTH_BY_FUNCTION_AND_STRESS_FIELD(
            RESU = RESU,
            INST = INST,
            GROUP_MA = GROUP_MA,
            GROUP_MA_CONT = GROUP_MA_CONT,
            MODELE = MODELE,
            MODELISATION = MODELISATION,
            MESH = MESH,
            GFUNC = GFUNC,
            FUNC_PARAMS = FUNC_PARAMS,
            SMOOTH_PARAMS = SMOOTH_PARAMS,
            GEOMETRIE = GEOMETRIE,
            SIGMA = __sigma)
    # }}}
    return ier
# Macro heading {{{
MORPHOGENESIS_GROWTH_BENEATH_CONTOUR = MACRO(
        nom = 'MORPHOGENESIS_GROWTH_BENEATH_CONTOUR',
        op = morphogenesis_growth_beneath_contour,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute the morphogenesis-inspired growth stress tensor.",
        RESU = SIMP(statut = "o", typ = resultat_sdaster),
        INST = SIMP(statut = "f", typ = "R"),
        GROUP_MA = SIMP(statut='o', typ = grma, validators = NoRepeat(),
            max='**'),
        GROUP_MA_CONT = SIMP(statut='o', typ = grma, validators = NoRepeat(),
            max='**'),
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        MODELISATION = SIMP(statut = "o", typ = 'TXM'),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        GFUNC = SIMP(statut = "d", typ = "TXM", defaut = 'mdf.Sgrowth'),
        FUNC_PARAMS = SIMP(statut = "f", typ = dict, max = "**"),
        SMOOTH_PARAMS = SIMP(statut = "o", typ = dict, max = "**"),
        GEOMETRIE = SIMP(statut = "d", typ = "TXM",
            defaut = 'DEFORMEE', into = ('DEFORMEE', 'INITIALE')),
        )
# }}}
# }}}

# Morphogenesis growth by inputting function and stress field {{{
def mg_growth_by_function_and_stress_field(self, RESU, INST, GROUP_MA,
        GROUP_MA_CONT, MODELE, MODELISATION, MESH, GFUNC, FUNC_PARAMS,
        SMOOTH_PARAMS, GEOMETRIE, SIGMA):
    # Initialisation {{{
    # Standard initialisation
    ier = 0
    self.set_icmd(1)
    # Define output variable
    self.DeclareOut("growth_t", self.sd)
    # Add required code_aster functions
    _F = self.get_cmd("_F")
    CREA_CHAMP = self.get_cmd("CREA_CHAMP")
    FORMULE = self.get_cmd("FORMULE")
    MODI_MAILLAGE = self.get_cmd("MODI_MAILLAGE")
    # }}}
    # Macro body {{{
    # Compute hydrostatic and shear stress {{{
    # Formulas {{{
    if MODELISATION == "D_PLAN":
        __Shyd_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ'),
                         VALE = 'Shyd(SIXX, SIYY, SIZZ)',
                         Shyd = mdf.Shyd_pln_strain)
        __Svon_f = FORMULE(NOM_PARA = ('SIXX', 'SIYY', 'SIZZ', 'SIXY'),
                         VALE = 'Svon(SIXX, SIYY, SIZZ, SIXY)',
                         Svon = mdf.Svmis_pln_strain)
    else:
        message  = "MODELISATION = '" + MODELISATION + "' has not been "
        message += "defined yet."
        raise ValueError(message)
    # }}}
    # Set up hydrostatic and shear stress as field {{{
    __hyshFu = CREA_CHAMP(OPERATION = "AFFE",
                          TYPE_CHAM = "NOEU_NEUT_F",
                          MODELE = MODELE,
                          AFFE = _F(TOUT = "OUI",
                                    NOM_CMP = ("X1", "X2"),
                                    VALE_F = (__Shyd_f, __Svon_f)))
    __hyshFi = CREA_CHAMP(OPERATION = "EVAL",
                          TYPE_CHAM = "NOEU_NEUT_R",
                          CHAM_F = __hyshFu,
                          CHAM_PARA = SIGMA)
    # }}}
    # }}}
    # Compute filtered-distance field {{{
    # Deform geometry if necessary {{{
    if GEOMETRIE == "DEFORMEE":
        # Get displacement field from results
        __depl = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                            OPERATION = 'EXTR',
                            RESULTAT = RESU,
                            NOM_CHAM = 'DEPL',
                            INST = INST)
        # Compute opposite displacement field
        __mdepl = CREA_CHAMP(TYPE_CHAM = 'NOEU_DEPL_R',
                             OPERATION = 'COMB',
                             COMB = _F(CHAM_GD = __depl,
                                    COEF_R = -1.0))
        # Deform the mesh
        MODI_MAILLAGE(reuse = MESH,
                      MAILLAGE = MESH,
                      DEFORME = _F(OPTION = 'TRAN',
                                   DEPL = __depl,),);
    # }}}
    # Set up the kdtree {{{
    __mail_py = MAIL_PY()
    __mail_py.FromAster(MESH)
    mail_coords = __mail_py.cn[:, :2]
    if isinstance(GROUP_MA_CONT, str):
         cont_nods = __mail_py.gno[GROUP_MA_CONT]
    else:
        cont_nods = []
        for group_ma_name in GROUP_MA_CONT:
            for nod_id in __mail_py.gno[group_ma_name]:
                cont_nods.append(nod_id)
    cont_coords = []
    for nod_id in cont_nods:
        cont_coords.append([mail_coords[nod_id, 0],
                            mail_coords[nod_id, 1]])
    cont_coords = np.array(cont_coords)
    cont_kdtree = cKDTree(cont_coords)
    # }}}
    # Distance-filter function {{{
    def DistanceFilter(xi, yi, kdtree, magnitude = 1.0):
        disti, _ = kdtree.query([xi, yi])
        smoothArg = SMOOTH_PARAMS[0]["smoothDis"]*(disti -
                SMOOTH_PARAMS[0]["maxDis"])
        if smoothArg > 705.0:
            return 0.0
        else:
            return magnitude/(1.0 + np.exp(smoothArg))
    DistanceFilter_cont = lambda xi, yi : \
            DistanceFilter(xi, yi, cont_kdtree, magnitude = 1.0)
    # }}}
    # Create code aster filtered-distance field {{{
    # Get geometry {{{
    __geom = CREA_CHAMP(TYPE_CHAM = 'NOEU_GEOM_R',
                        NOM_CHAM = 'GEOMETRIE',
                        OPERATION = "EXTR",
                        MAILLAGE = MESH)
    if GEOMETRIE == "DEFORMEE":
        # Deform the mesh backwards
        MODI_MAILLAGE(reuse = MESH,
                      MAILLAGE = MESH,
                      DEFORME = _F(OPTION = 'TRAN',
                                   DEPL = __mdepl,),);
    # }}}
    __distFo = FORMULE(NOM_PARA = ("X", "Y"),
                       VALE = "FUNC(X, Y)",
                       FUNC = DistanceFilter_cont)
    __distFu = CREA_CHAMP(OPERATION = "AFFE",
                          TYPE_CHAM = "NOEU_NEUT_F",
                          MODELE = MODELE,
                          AFFE = _F(TOUT = "OUI",
                                    NOM_CMP = "X3",
                                    VALE_F = __distFo))
    __distFi = CREA_CHAMP(OPERATION = "EVAL",
                          TYPE_CHAM = "NOEU_NEUT_R",
                          CHAM_F = __distFu,
                          CHAM_PARA = __geom)
    # }}}
    # }}}
    # Compute filtered stresses {{{
    __filFoH = FORMULE(NOM_PARA = ("X1", "X3"),
                       VALE = "X1*X3")
    __filFoS = FORMULE(NOM_PARA = ("X2", "X3"),
                       VALE = "X2*X3")
    __filFun = CREA_CHAMP(OPERATION = "AFFE",
                          TYPE_CHAM = "NOEU_NEUT_F",
                          MODELE = MODELE,
                          AFFE = _F(TOUT = "OUI",
                                    NOM_CMP = ("X1", "X2"),
                                    VALE_F = (__filFoH, __filFoS)))
    __filStr = CREA_CHAMP(OPERATION = "EVAL",
                          TYPE_CHAM = "NOEU_NEUT_R",
                          CHAM_F = __filFun,
                          CHAM_PARA = [__hyshFi, __distFi])
    # }}}
    # Compute growth stress tensor {{{
    # Maximum stresses
    tauMiswp = __filStr.EXTR_COMP('X2', []).valeurs
    sigHydwp = __filStr.EXTR_COMP('X1', []).valeurs
    maxSigHydCom = abs(min(sigHydwp))
    maxTau = max(tauMiswp)
    # Formula
    sgfunc = eval(GFUNC)
    func = lambda X1, X2: sgfunc(X1/maxSigHydCom, X2/maxTau,
            **FUNC_PARAMS[0])
    __Sgrw_f  = FORMULE(NOM_PARA = ('X1', 'X2'),
                        VALE = 'func(X1, X2)',
                        func = func,
                        )
    # Function field
    __grwFu = CREA_CHAMP(OPERATION = "AFFE",
                          TYPE_CHAM = "NOEU_NEUT_F",
                          MODELE = MODELE,
                          AFFE = _F(GROUP_MA = GROUP_MA,
                                    NOM_CMP = ("X1"),
                                    VALE_F = (__Sgrw_f)))
    # Field
    __grw_fi = CREA_CHAMP(OPERATION = "EVAL",
                          TYPE_CHAM = "NOEU_NEUT_R",
                          CHAM_F = __grwFu,
                          CHAM_PARA = __filStr)
    # Growth stress as tensor
    growth_t = CREA_CHAMP(OPERATION = 'ASSE',
                          TYPE_CHAM = 'NOEU_SIEF_R',
                          MODELE = MODELE,
                          ASSE = (_F(GROUP_MA = GROUP_MA,
                                    CHAM_GD = __grw_fi,
                                    NOM_CMP = ('X1'),
                                    NOM_CMP_RESU = ('SIXX')),
                                  _F(GROUP_MA = GROUP_MA,
                                    CHAM_GD = __grw_fi,
                                    NOM_CMP = ('X1'),
                                    NOM_CMP_RESU = ('SIYY')),
                                  _F(GROUP_MA = GROUP_MA,
                                    CHAM_GD = __grw_fi,
                                    NOM_CMP = ('X1'),
                                    NOM_CMP_RESU = ('SIZZ')),),)
    # }}}
    # }}}
    return ier
# Macro heading {{{
MG_GROWTH_BY_FUNCTION_AND_STRESS_FIELD = MACRO(
        nom = 'MG_GROWTH_BY_FUNCTION_AND_STRESS_FIELD',
        op = mg_growth_by_function_and_stress_field,
        sd_prod = cham_no_sdaster,
        docu = "Macro to compute the morphogenesis-inspired growth stress tensor.",
        RESU = SIMP(statut = "o", typ = resultat_sdaster),
        SIGMA = SIMP(statut = "o", typ = cham_no_sdaster),
        INST = SIMP(statut = "f", typ = "R"),
        GROUP_MA = SIMP(statut='o', typ = grma, validators = NoRepeat(),
            max='**'),
        GROUP_MA_CONT = SIMP(statut='o', typ = grma, validators = NoRepeat(),
            max='**'),
        MODELE = SIMP(statut='o', typ=modele_sdaster),
        MODELISATION = SIMP(statut = "o", typ = 'TXM'),
        MESH = SIMP(statut = "o", typ = maillage_sdaster),
        GFUNC = SIMP(statut = "d", typ = "TXM", defaut = 'mdf.Sgrowth'),
        FUNC_PARAMS = SIMP(statut = "f", typ = dict, max = "**"),
        SMOOTH_PARAMS = SIMP(statut = "o", typ = dict, max = "**"),
        GEOMETRIE = SIMP(statut = "d", typ = "TXM",
            defaut = 'DEFORMEE', into = ('DEFORMEE', 'INITIALE')),
        )
# }}}
# }}}
# }}}

# Functions {{{
# Save contours {{{
def SaveContour(Unit, fileName, groupName, inst = None):
    from Noyau.N__F import _F
    from code_aster.Cata.Commands import DEFI_FICHIER, IMPR_RESU
    DEFI_FICHIER(ACTION = 'ASSOCIER',
                 FICHIER = fileName,
                 UNITE = Unit)
    if inst == None:
        IMPR_RESU(RESU=_F(RESULTAT = resu,
                          NOM_CHAM = 'DEPL',
                          #INST = 1,
                          GROUP_MA = groupName),
                  FORMAT = 'GMSH',
                  UNITE=Unit)
    else:
        IMPR_RESU(RESU=_F(RESULTAT = resu,
                          NOM_CHAM = 'DEPL',
                          INST = inst,
                          GROUP_MA = groupName),
                  FORMAT = 'GMSH',
                  UNITE=Unit)
    DEFI_FICHIER(ACTION = 'LIBERER', UNITE = U)
    # }}}

# Nodal wear displacements _F {{{
def NodalWearDisplacements_F(wear_status, mail_py):
    from Noyau.N__F import _F
    FFF = []
    for k1 in range(len(wear_status)):
        nod_i = mail_py.correspondance_noeuds[wear_status[k1][0] - 1]
        ux_i = wear_status[k1][1]
        uy_i = wear_status[k1][2]
        FFF.append(_F(NOEUD = nod_i,
                      DX = ux_i,
                      DY = uy_i))
    return FFF
# }}}

# Element equivalent pressure {{{
def ElementEquivalentPressure_F(x, p_c, co_info, mail_py, area = -1.0):
    from Noyau.N__F import _F
    # Active elements
    activeElements = co_info.EXTR_TABLE().values()["ACTIVE_ELEMENTS"]
    try:
        activeElements.remove(None)
    except ValueError:
        pass
    # Contact area
    if area <= 0.0:
        contArea = co_info.EXTR_TABLE().values()["CONT_AREA"][0]
    else:
        contArea = area
        firstElement = co_info.EXTR_TABLE().values()["OUT_ELEMENTS"][0]
        activeElements = mdf.OrganiseElementSet(mail_py, activeElements, firstElement)
        nodes = mail_py.cn[:, :2]
        eles = mail_py.co
    p_o = abs(simps(p_c, x)/contArea)
    # Active element names
    activeElementNames = []
    accarea = 0.0
    for idd in activeElements:
        if area > 0.0:
            ele = eles[idd - 1]
            x1 = nodes[ele[0], 0]
            x2 = nodes[ele[1], 0]
            y1 = nodes[ele[0], 1]
            y2 = nodes[ele[1], 1]
            accarea += np.sqrt((x2 - x1)**2.0 + (y2 - y1)**2.0)
            if accarea <= contArea:
                activeElementNames.append(mail_py.correspondance_mailles[idd - 1])
        else:
            activeElementNames.append(mail_py.correspondance_mailles[idd - 1])
    # Assign pressure to active elements
    FFF = []
    for eleName in activeElementNames:
        FFF.append(_F(PRES = p_o,
                      MAILLE = eleName))
    return FFF
# }}}
# }}}
