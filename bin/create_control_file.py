namelist_control = [ \
                     'lnmlj',\
                     'lcoulomb',\
                     'lmorse',\
                     'lbmhft',\
                     'lbmhftd',\
                     'lsurf',\
                     'lcsvr',\
                     'lharm',\
                     'ltraj',\
                     'lvnlist',\
                     'lstatic',\
                     'lreduced',\
                     'lreducedN',\
                     'ltest',\
                     'lmsd',\
                     'lvacf',\
                     'lrestart',\
                     'cutlongrange',\
                     'cutshortrange',\
                     'calc',\
                     'dgauss',\
                     'longrange',\
                     'itraj_start',\
                     'itraj_period',\
                     'itraj_format',\
                     'trajff_data',\
                     'iscff_format',\
                     'iscff_data',\
                     'iefall_format',\
                     'iefgall_format',\
                     'idipall_format',\
                     'restart_data',\
                     'skindiff']
namelist_md = [ \
                     'integrator',\
                     'setvel',\
                     'npas',\
                     'nequil',\
                     'nequil_period',\
                     'annealing',\
                     'npropr',\
                     'npropr_start',\
                     'nprint',\
                     'nprint',\
                     'fprint',\
                     'spas',\
                     'dt ',\
                     'temp',\
                     'press',\
                     'nuandersen',\
                     'taucsvr',\
                     'tauTberendsen',\
                     'tauPberendsen',\
                     'nhc_yosh_order',\
                     'nhc_mults',\
                     'nhc_n',\
                     'timesca_thermo',\
                     'timesca_baro' ]

namelist_field = [ \
                     'lKA',\
                     'ctrunc',\
                     'symmetric_pot',\
                     'ncelldirect',\
                     'kES',\
                     'alphaES',\
                     'qlj',\
                     'plj',\
                     'sigmalj',\
                     'epslj',\
                     'sigmamor',\
                     'epsmor',\
                     'rhomor',\
                     'mass',\
                     'doefield',\
                     'doefg',\
                     'qch',\
                     'quad_nuc',\
                     'dip',\
                     'quad',\
                     'poldip',\
                     'poldip_iso',\
                     'polquad',\
                     'polquad_iso',\
                     'pol_damp_b',\
                     'pol_damp_c',\
                     'pol_damp_k',\
                     'extrapolate_order',\
                     'conv_tol_ind',\
                     'min_scf_pol_iter',\
                     'max_scf_pol_iter',\
                     'algo_moment_from_pola',\
                     'algo_ext_dipole',\
                     'thole_functions',\
                     'thole_function_type',\
                     'thole_param',\
                     'omegakO',\
                     'epsw',\
                     'lautoES',\
                     'lwfc',\
                     'lwrite_dip_wfc',\
                     'lwrite_dip',\
                     'lwrite_quad',\
                     'lwrite_ef',\
                     'lwrite_efg',\
                     'ldip_wfc',\
                     'rcut_wfc',\
                     'ldip_polar',\
                     'ldip_damping',\
                     'lquad_polar',\
                     'Abmhftd',\
                     'Bbmhftd',\
                     'Cbmhftd',\
                     'Dbmhftd',\
                     'BDbmhftd']

namelist_stochio = [ 'def',\
                     'typedef',\
                     'density',\
                     'lcubic',\
                     'a_o_b',\
                     'a_o_c',\
                     'target_nions',\
                     'atoms_in',\
                     'natoms_in',\
                     'sio2',\
                     'geo2',\
                     'na2o',\
                     'b2o3',\
                     'la2o3',\
                     'sro',\
                     'cao',\
                     'al2o3',\
                     'moo3',\
                     'p2o5'] 

namelist_all = namelist_md + namelist_field + namelist_control + namelist_stochio

calc_allowed = [ "'md'" , "'opt'" , "'vib'" , "'vib+fvib'" , "'vib+gmod'" , \
                 "'vib+band'" , "'vib+dos'" , "'efg'" , "'efg+acf'" , "'gr'" ,  "'vois1'" , "'rmc'" , "'dist'" ,"'stochio'"]
longrange_allowed = [ "'ewald'" , "'direct'" ]
dgauss_allowed = [ "'boxmuller_basic'", "'boxmuller_polar'" , "'knuth'" ]
data_allowed = [ "'rvf'" , "'rnn'" , "'rnf'" , "'rvn'" ]
ctrunc_allowed = [ "'notrunc'", "'linear'" , "'quadratic'" ]
thole_function_type_allowed = [ "'linear'","'expon1'","'expon2'", "'gauss'" ]
algo_ext_dipole_allowed  = [ "'poly'","'aspc'" ]
algo_moment_from_pola_allowed =[ "'scf'" , "'scf_kO_v1'" , "'scf_kO_v2'" , "'scf_kO_v3'" , "'scf_kO_v4_1'" , "'scf_kO_v4_2'"]
integrator_allowed = [ "'nve-vv'"  , "'nve-lf'"   , "'nve-be'" ,  "'nve-lfq'", "'nvt-and'" , "'nvt-nhc2'" , "'nvt-nhcn'", "'npe-vv'"  , "'npt-nhcnp'" ]
setvel_allowed = [ "'MaxwBoltz'", "'Uniform'" ]
typedef_allowed = [ 'oxydes' , 'atomic' ]
def_allowed = [ 'pct' , 'num' ]

all_values_allowed = calc_allowed + longrange_allowed + dgauss_allowed + data_allowed + \
                     ctrunc_allowed + thole_function_type_allowed + algo_ext_dipole_allowed + \
                     algo_moment_from_pola_allowed + integrator_allowed + setvel_allowed + typedef_allowed + def_allowed

default_values={}

default_values['lnmlj'] = ".false."
default_values['lbmhft'] = ".false."
default_values['lbmhftd'] = ".false."
default_values['lmorse'] = ".false."
default_values['lcoulomb'] = ".false."
default_values['lsurf'] = ".false."
default_values['lcsvr'] = ".false."
default_values['lharm'] = ".false."
default_values['ltraj'] = ".false."
default_values['lvnlist'] = ".true."
default_values['lstatic'] = ".false."
default_values['lreduced'] = ".false."
default_values['lreducedN'] = ".false."
default_values['ltest'] = ".false."
default_values['lmsd'] = ".false."
default_values['lvacf'] = ".false."
default_values['lrestart'] = ".false."
default_values['full_restart'] = ".false."
default_values['calc'] = "'md'"
default_values['dgauss'] = "'boxmuller_basic'"
default_values['longrange'] = "'ewald'"
default_values['skindiff'] = "0.15_dp"
default_values['cutshortrange'] = "0.0_dp"
default_values['cutlongrange'] = "0.0_dp"
default_values['itraj_start'] = "1"
default_values['itraj_period'] = "10000"
default_values['itraj_format'] = "1"
default_values['iscff_format'] = "1"
default_values['iefall_format'] = "1"
default_values['iefgall_format'] = "1"
default_values['idipall_format'] = "1"
default_values['trajff_data'] = "'rnn'"
default_values['iscff_data'] = "'rnn'"
default_values['restart_data'] = "'rnn'"
default_values['itime'] = "0"
default_values['lleapequi'] = ".false."
default_values['integrator'] = "'nve-vv'"
default_values['setvel'] = "'MaxwBoltz'"
default_values['npas'] = "10"
default_values['nequil'] = "0"
default_values['nequil_period'] = "1"
default_values['nprint'] = "1"
default_values['fprint'] = "1"
default_values['spas'] = "1000"
default_values['dt '] = "0.0_dp"
default_values['lKA'] = ".false."
default_values['temp'] = "1.0_dp"
default_values['press'] = "0.0_dp"
default_values['nuandersen'] = "0.0_dp"
default_values['taucsvr'] = "0.0_dp"
default_values['tauTberendsen'] = "0.0_dp"
default_values['tauPberendsen'] = "0.0_dp"
default_values['nhc_yosh_order'] = "3"
default_values['nhc_mults'] = "2"
default_values['nhc_n'] = "4"
default_values['annealing'] = "1.0_dp"
default_values['npropr'] = "1"
default_values['npropr_start'] = "0"
default_values['timesca_thermo'] = "1.0_dp"
default_values['timesca_baro'] = "1.0_dp"
default_values['ctrunc'] = "'notrunc'"
default_values['symmetric_pot'] = ".true."
default_values['lKA'] = ".false."
default_values['epslj'] = "0.0_dp"
default_values['sigmalj'] = "0.0_dp"
default_values['qlj'] = "12.0_dp"
default_values['plj'] = "6.0_dp"
default_values['epsmor'] = "0.0_dp"
default_values['sigmamor'] = "0.0_dp"
default_values['rhomor'] = "0.0_dp"
default_values['Abmhftd'] = "0.0_dp"
default_values['Bbmhftd'] = "0.0_dp"
default_values['Cbmhftd'] = "0.0_dp"
default_values['Dbmhftd'] = "0.0_dp"
default_values['BDbmhftd'] = "0.0_dp"
default_values['ncelldirect'] = "2"
default_values['kES'] = "10 10 10"
default_values['alphaES'] = "1.0_dp"
default_values['epsw'] = "1e-6"
default_values['lautoES'] = ".false."
default_values['qch'] = "0.0_dp"
default_values['quad_nuc'] = "0.0_dp"
default_values['dip'] = "0.0_dp"
default_values['quad'] = "0.0_dp"
default_values['doefield'] = ".false."
default_values['doefg'] = ".false."
default_values['lwrite_dip'] = ".false."
default_values['lwrite_quad'] = ".false."
default_values['lwrite_ef'] = ".false."
default_values['lwrite_efg'] = ".false."
default_values['ldip_polar'] = ".false."
default_values['lquad_polar'] = ".false."
default_values['poldip'] = "0.0_dp"
default_values['poldip_iso'] = "0.0_dp"
default_values['polquad'] = "0.0_dp"
default_values['polquad_iso'] = "0.0_dp"
default_values['pol_damp_b'] = "0.0_dp"
default_values['pol_damp_c'] = "0.0_dp"
default_values['pol_damp_k'] = "0"
default_values['ldip_damping'] = ".false."
default_values['conv_tol_ind'] = "1e-6"
default_values['min_scf_pol_iter'] = "3"
default_values['max_scf_pol_iter'] = "100"
default_values['extrapolate_order'] = "0"
default_values['algo_ext_dipole'] = "'aspc'"
default_values['algo_moment_from_pola'] = "'scf'"
default_values['thole_functions'] = ".false."
default_values['thole_function_type'] = "'linear'"
default_values['thole_param'] = "1.662_dp"
default_values['omegakO'] = "0.7_dp"
default_values['lwfc'] = "0"
default_values['lwrite_dip_wfc'] = ".false."
default_values['ldip_wfc'] = ".true."
default_values['rcut_wfc'] = "0.5_dp"
default_values['mass'] = "1.0_dp"
default_values['task_coul'] = ".false."
default_values['typedef'] = "'oxydes'"
default_values['lcubic'] = ".true."
default_values['a_o_b'] = "1.0_dp"
default_values['a_o_c'] = "1.0_dp"
default_values['sio2'] = "0._dp"
default_values['na2o'] = "0._dp"
default_values['b2o3'] = "0._dp"
default_values['la2o3'] = "0._dp"
default_values['sro'] = "0._dp"
default_values['cao'] = "0._dp"
default_values['p2o5'] = "0._dp"
default_values['al2o3'] = "0._dp"
default_values['moo3'] = "0._dp"
default_values['geo2'] = "0._dp"
default_values['def'] = "'pct'"
default_values['density'] = "0.0d0"
default_values['atoms_in'] = ""
default_values['natoms_in'] = "0"
default_values['target_nions'] = "0"

def get_control_file(calc,md_values):

    tags=['control','md','field','stochio']

    tag_namelist={}
    tag_namelist['control'] = namelist_control
    tag_namelist['md'] = namelist_md
    tag_namelist['field'] = namelist_field
    tag_namelist['stochio'] = namelist_stochio

    tags_md=['control','md','field']


    filename = "control_"+calc+".F"
    file = open(filename,"w")
    file.write("!"+60*"="+"\n")
    file.write("!   Control file generated by MDff GUI\n")
    file.write("!   generator key: "+calc+"\n")
    file.write("!"+60*"="+"\n")

    if calc == "full" :

        for tag in tags:
            file.write("&"+tag+"tag\n")
            file.write("")
            for namelist in tag_namelist[tag]:
                file.write("        "+namelist+" = "+default_values[namelist]+"\n")
            file.write("&end\n")
            file.write("\n")
                
    if calc == "empty" :

        for tag in tags:
            file.write("&"+tag+"tag\n")
            file.write("&end\n")
            file.write("\n")

    if calc == "stochio" :
       
        file.write("&controltag\n")
        file.write("        calc = 'stochio'\n")
        file.write("&end\n")
        file.write("")
        file.write("&stochiotag\n")
        file.write("")
        for namelist in tag_namelist['stochio']:
            file.write("        "+namelist+" = "+default_values[namelist]+"\n")
        file.write("&end\n")
        file.write("\n")

    if calc == "md":

        file.write("&controltag\n")
        file.write("        calc = 'md'\n")
        file.write("&end\n")
        file.write("")
        file.write("&md\n")
        file.write("")
        for namelist in tag_namelist['md']:
            file.write("        "+namelist+" = "+md_values[namelist]+"\n")
        file.write("&end\n")
        file.write("\n")

    file.close()

        



if __name__ == "__main__":

    get_control_file("full")
