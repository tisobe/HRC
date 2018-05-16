#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################################
#                                                                                                                           #
#   OcatSQL.py: read data from sql database to fill all Ocat Data Page pamater values                                       #
#                                                                                                                           #
#       author: t. isobe (tisobe@cfa.harvard.edu)                                                                           #
#                                                                                                                           #
#       last update: Jan 06, 2017                                                                                           #
#                                                                                                                           #
#############################################################################################################################

import sys
import os
import string
import re
import copy
import unittest
#
#--- reading directory list
#
path = '/data/mta4/CUS/www/Usint/Ocat_sub/house_keeping/dir_list_py'

f    = open(path, 'r')
data = [line.strip() for line in f.readlines()]
f.close()

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec "%s = %s" %(var, line)
#
#--- append path to a private folders
#
sys.path.append(bin_dir)
#sys.path.append(mta_dir)
#
#--- sybase module
#
import DBI as dbi

#----------------------------------------------------------------------------------------------------------------
#-- OcatDB: create an Ocat Data Object                                                                        ---
#----------------------------------------------------------------------------------------------------------------

class OcatDB(object):
    """
    create an Ocat Data Object for a given obsid
    Usage:
    db = OcatDB(15113)

    db.paramList()                       --- list all parameter names
    db.origValue('y_det_offset')         --- give the original value of the parameterr ('y_det_offset ')
                                             if option ordr=<order #> is given, show the value for the order
    db.newValue('y_det_offset')          --- give the current value of the parameter
                                             if option ordr=<order #> is given, show the value for the order
    db.updateValue('y_det_offset',0.14)  --- update the parameter value
                                             if option ordr=<order #> is given, update the value for the order
    db.addParam('newparam' value='1234') --- add a new parameter with an optional value. without option 
                                             the value is ""
    db.addNewOrder('roll_ordr', ['P', 'Y', '90', '14])
                                         --- add new order to the given order possible parameters
                                             the value list must follow extact order of the paramters
                                             see <hosue_keeping>/ocat_name_table for details 
                                             (roll_rodr: roll_table/ time_rodr: time_table/ ordr: aciswin_table)
    db.getTableEntries('acis_table')    ---  give back parameter names for a given table group

    """
#---------------------------------------------------------

    def __init__(self, obsid):
        self.db = get_ocat_param(str(obsid))

#---------------------------------------------------------

    def paramList(self):
        """
        return all parameter names in the database
        """
        return self.db.keys()

#---------------------------------------------------------

    def origValue(self, param, ordr=''):
        """
        return the original value of the parameter
        """
        try:
            if ordr == '':
                return self.db[str(param)][0] 
            else:
                rank = int(float(ordr)) -1
                return self.db[str(param)][0][rank] 
        except:
            return 'NA'

#---------------------------------------------------------

    def newValue(self, param, ordr=''):
        """
        return the current value of the parameter
        """
        try:
            if ordr == '':
                return self.db[str(param)][1] 
            else:
                rank = int(float(ordr)) -1
                return self.db[str(param)][1][rank]
        except:
            return 'NA'

#---------------------------------------------------------

    def updateValue(self, param, val, ordr=''):
        """
        update the current value. if ordr is given, the <ordr>_th value of the parameter will be updated
        """
        if ordr == '':
            self.db[str(param)][1] = val 
        else:
            try:
                pos = int(float(ordr)) -1
                if pos < 0:
                    pos = 0
            except:
                pos = 0
            try:
                self.db[str(param)][1][pos] = val           #--- for the case, exchanging the value
            except:
                self.db[str(param)][1].append(val)          #--- for the case, adding a new value

#---------------------------------------------------------

    def addParam(self, param, value=''):
        """
        add parameter and its value to the object database
        """
        self.db[param]    = []
        if value != '':
            self.db[param] = [value, value]

#---------------------------------------------------------

    def addNewOrder(self, oname, newlist):
        """
        add new order to the given parameter group
        Input: oname: roll_odr, time_ordr, ordr (aciswin)
               newlist: a list of values. Values must be given in the following order
               roll_table    : roll_constraint, roll_180, roll, roll_tolerance
               time_table    : window_constraint, tstart, tstop
               aciswin_table : start_row, start_column, width, height, lower_threshold, pha_range, sample, chip, include_flag , aciswin_id
        """

        table_dict = read_db_param_table()
        if oname == 'roll_ordr':
            plist = table_dict['roll_table'][1]
        elif oname == 'time_ordr':
            plist = table_dict['time_table'][1]
        else:
            plist = table_dict['aciswin_table'][1]

        try:
            ordr = int(float(self.newValue(oname)) + 1)
        except:
            ordr = 1

        self.updateValue(oname, ordr)


        for  i in range(0, len(plist)):
            try:
                vlist  = self.newValue(plist[i])
            except:
                vlist = []
            vlist.append(newlist[i])

            self.updateValue(plist[i], vlist)

#---------------------------------------------------------

    def removeNullOrder(self, oname):
        """
        remove "NULL" rows from the ordered entries. 
        Input:  oname: roll_odr, time_ordr, ordr (aciswin)
        """

        table_dict = read_db_param_table()
        if oname == 'roll_ordr':
            plist = table_dict['roll_table'][1]
            lead  = 'roll_constraint'
            flag  = 'roll_flag'
        elif oname == 'time_ordr':
            plist = table_dict['time_table'][1]
            lead  = 'window_constraint'
            flag  = 'window_flag'
        else:
            plist = table_dict['aciswin_table'][1]
            lead  = 'chip'
            flag  = 'spwindow_flag'
#
#--- check what is the current order. If it is not defined, set to 1
#
        try:
            ordr = int(float(self.newValue(oname)))
        except:
            ordr = 1
#
#--- find which order is NULL
#
        remove_list = []
        for i in range(0, ordr):
            rank = i + 1
            val = self.newValue(lead, ordr=rank)
            if val == 'NULL' or val == 'NO':
                remove_list.append(i)
#
#--- removed the entries of i-th order
#
        offset = 0                                  #---- monitor how many rows are removed
        for i in remove_list:
            k = i - offset
            offset += 1
            for ent in plist:
                vlist = self.newValue(ent)
                del vlist[k]
                self.updateValue(ent, vlist)
#
#--- update the order value
#        
        ordr -= offset
        if ordr < 0:
            ordr = 0

        self.updateValue(oname, ordr)
#
#--- if there is no more order entry, set flag to "NULL"
#
        if ordr == 0:
            self.updateValue(flag, "NULL")

#---------------------------------------------------------

    def getTableEntries(self, table_name):
        """
        give back parameter names for a given table group
        Input:  table_name: e.g. "aciswin_table"
        """

        table_dict = read_db_param_table()
        return table_dict[table_name][1]

#----------------------------------------------------------------------------------------------------------------
#--- get_ocat_param: extract a basic target information for a given obsid                                    ---
#----------------------------------------------------------------------------------------------------------------

def get_ocat_param(obsid):

    """
    extract all parameters needed for Ocat Data Page for a given obsid.

    Input:  obsid       obsid of the observation
    Output: db_dict     dictionary of extracted parameters and their values
                        db_dict ={<parameter name>, [<value from DB>,<updatable value>]}

    Note: we aslo need to access <house_keeping>/ocat_name_list_table
    """
#
#--- read lists of parameters name
#
    table_dict = read_db_param_table()
#
#--- start filling parmaters/ values in the dictionary form 
#
    db_dict = {'obsid' : [obsid, obsid]}
#
#--- read general information (target_table) parameters
#
    criteria = ' from target where obsid=' + str(obsid)
    temp_dict = extract('axafocat', table_dict, criteria,'target_table')
    db_dict.update(temp_dict)
#
#--- read triggering related flag information (flag_table) parameters
#
    criteria = ' from target where obsid=' + str(obsid)
    temp_dict = extract('axafocat', table_dict, criteria,'flag_table')
    db_dict.update(temp_dict)
#
#--- read dither (dither_table) parameters
#
    criteria = ' from dither where  obsid=' + str(obsid)
    temp_dict = extract('axafocat', table_dict, criteria, 'dither_table')
    db_dict.update(temp_dict)
#
#--- read sim (sim_table) paremeters
#
    criteria = ' from sim where  obsid=' + str(obsid)
    temp_dict = extract('axafocat', table_dict, criteria, 'sim_table')
    db_dict.update(temp_dict)
#
#--- read soe (soe_table) parameter
#
    criteria = ' from soe where  obsid=' + str(obsid)
    temp_dict = extract('axafocat', table_dict, criteria, 'soe_table')
    db_dict.update(temp_dict)
#
#--- read proposal related (prop_table) parameters
#
    criteria = ' from prop_info where  ocat_propid=' + str(db_dict['ocat_propid'][0])
    temp_dict = extract('axafocat', table_dict, criteria, 'prop_table')
    db_dict.update(temp_dict)
#
#--- roll constraints (roll_table) parameters
#
    criteria   = 'from rollreq where obsid=' + str(obsid) + ' and  ordr='
    temp_dict = extract_order_entry('axafocat', table_dict, criteria, 'roll_table', 'roll_ordr')
    db_dict.update(temp_dict)
    try:
        val = db_dict['roll_ordr']
    except:
        db_dict['roll_ordr'] = [0, 0]
#
#--- time constraints (time_table) parameters
#
    criteria   = 'from timereq where obsid=' + str(obsid) + '  and ordr='
    temp_dict = extract_order_entry('axafocat', table_dict, criteria, 'time_table', 'time_ordr')
    db_dict.update(temp_dict)
    try:
        val = db_dict['time_ordr']
    except:
        db_dict['time_ordr'] = [0, 0]
#
#--- read other costratint (const_table) parameters
#
    criteria = ' from target where  obsid=' + str(obsid)
    temp_dict = extract('axafocat', table_dict, criteria, 'const_table')
    db_dict.update(temp_dict)
#
#--- read phase costratint (phase_table) parameters
#
    criteria = ' from phasereq where  obsid=' + str(obsid)
    temp_dict = extract('axafocat', table_dict, criteria, 'phase_table')
    db_dict.update(temp_dict)
#
#--- read too costratint (too_table) parameters
#
    if db_dict['tooid'][0]:
        criteria = ' from too where  tooid=' + str(db_dict['tooid'][0])
        temp_dict = extract('axafocat', table_dict, criteria, 'too_table')
        db_dict.update(temp_dict)
    else:
        for ent in table_dict['too_table'][1]:
            db_dict[ent] = ['', '']
#
#--- acis parameters (acis_table)
#
    if db_dict['acisid'][0]:
        criteria = ' from acisparam where  acisid=' + str(db_dict['acisid'][0])
        temp_dict = extract('axafocat', table_dict, criteria, 'acis_table')
        db_dict.update(temp_dict)
    else:
        for ent in table_dict['acis_table'][1]:
            db_dict[ent] = ['', '']
#
#--- acis window constraint (aciswin_id_table/aciswin_table) parameters 
#
    criteria = ' from aciswin where  obsid=' + str(obsid)
    temp_dict = extract('axafocat', table_dict, criteria, 'aciswin_id_table')
    db_dict.update(temp_dict)

    if db_dict['aciswin_id'][0]:
        criteria   = 'from aciswin where obsid = ' +  str(obsid) + ' and ordr='
        temp_dict = extract_order_entry('axafocat', table_dict, criteria, 'aciswin_table', 'ordr')
        db_dict.update(temp_dict)
        awc_l_th = 0
        try:
            if db_dict['ordr'][0] > 0:
                for ent in db_dict['lower_threshold'][0]:
                    if ent > 0.5:
                        awc_l_th  = 1
        except:
            pass
    
        db_dict['awc_l_th'] = [awc_l_th, awc_l_th]
    else:
        for ent in table_dict['aciswin_table'][1]:
#            db_dict[ent] = ['', '']
            db_dict[ent] = [[],[]]
    try:
        val = db_dict['ordr']
    except:
        db_dict['ordr'] = [0, 0]
#
#--- hrc parameters (hrc_table)
#
    if db_dict['hrcid'][0]:
        criteria = ' from hrcparam where  hrcid=' + str(db_dict['hrcid'][0])
        temp_dict = extract('axafocat', table_dict, criteria, 'hrc_table')
        db_dict.update(temp_dict)
    else:
        for ent in table_dict['hrc_table'][1]:
            db_dict[ent] = ['', '']
#
#--- checking monitoring status
#
    if db_dict['pre_id'][0]:
        monitor_flag = 'Y'
    else:
        monitor_flag = 'N'
        try:
            cmd  = 'select distinct pre_id from target where pre_id=' + str(obsid)
            out  = readSQL(cmd, 'axafocat')
    	    if out['pre_id'] is not None:
        	    monitor_flag = 'Y'
        except:
	        pass
#
#--- check group entries
#
    if db_dict['group_id'][0]:
        monitor_flag = 'N'
        db_dict['pre_min_lead'] = ['None', 'None']
        db_dict['pre_max_lead'] = ['None', 'None']
        db_dict['pre_id']       = ['None', 'None']
        db_dict['monitor']      = [[], []]

        cmd = 'select obsid from target where group_id="' + db_dict['group_id'][0] + '"'
        out = readSQL(cmd, 'axafocat', all='yes')
        group = []
        for ent in out:
            val = ent[0]
            group.append(val)

        db_dict['group'] = [group, group]
    else:
        db_dict['group'] = [[],[]]
#
#--- if monitor flag is still "Y", find which obsids blong to this monitor list
#
    if monitor_flag == 'Y':
        vent  = find_monitor_obs(obsid)
        db_dict['monitor']  = [vent, vent]
        db_dict['monitor_flag']  = ['Y', 'Y']
    else:
        db_dict['monitor']  = [[],[]]
        db_dict['monitor_flag']  = ['N', 'N']
#
#--- update ao #
#
    cmd = 'select ao_str from prop_info where ocat_propid=' + str(db_dict['ocat_propid'][0])
    out = readSQL(cmd, 'axafocat')
    db_dict['obs_ao_str'] = [out['ao_str'], out['ao_str']]
    db_dict.update(temp_dict)
#
#--- proporser's name
#
    criteria = ' from person_short s,axafocat..prop_info p where  p.ocat_propid=' + str(db_dict['ocat_propid'][0]) + ' and s.pers_id=p.piid'
    temp_dict = extract('axafusers', table_dict, criteria, 'pi_table')
    db_dict.update(temp_dict)
#
#--- check co-investigator (observer)
#
    criteria  = ' from person_short s,axafocat..prop_info p where  p.ocat_propid=' + str(db_dict['ocat_propid'][0]) + ' and s.pers_id=p.piid'
    temp_dict =  extract('axafusers', table_dict, criteria, 'co_i_table')

    mc = re.search('Y|y', temp_dict['observer'][0])
    if mc is not None:
        criteria = ' from person_short s,axafocat..prop_info p where  p.ocat_propid=' + str(db_dict['ocat_propid'][0]) + ' and p.coin_id = s.pers_id'
        temp_dict = extract('axafusers',  table_dict, criteria, 'observer_table')
        db_dict.update(temp_dict)
    else:
        db_dict['observer'] = db_dict['pi_name']
#
#--- read remarks (remark_table) 
#
    criteria = ' from target where obsid=' + str(obsid)
    temp_dict = extract('axafocat', table_dict, criteria,'remark_table')
    db_dict.update(temp_dict)
#
#--- add parameters which are not from database
#
    planned_roll            = find_planned_roll(str(obsid))
    db_dict['planned_roll'] = [planned_roll, planned_roll]
    db_dict['comments']     = ['', '']

    return db_dict

#----------------------------------------------------------------------------------------------------------------
#-- read_db_param_table: read ocat_name_table and create a dictionary of parameter lists                      ---
#----------------------------------------------------------------------------------------------------------------

def read_db_param_table():

    """
    read "ocat_name_table" and create a dictionary of parameter lists
        NOTE: if you need to add parameters to be read from sql database, add to <house_keeping>ocat_name_table

    Input:  <hosue_keeping>/ocat_name_table
    Output: table_dict in the form of {<table_name>:[<db parameter name list>,<ocat data parameter list>]
            if <ocat data parameter list> is empty at ocat_name_table, <db parameter name list> is used
    """
#
#--- ocat_name_table contains <table_name> and list of parameters related to the table
#
    file = house_keeping + '/ocat_name_table'
    f    = open(file, 'r')
    data = [line.strip() for line in f.readlines()]
    f.close()

    table_dict = {}
    for ent in data:
        if ent.startswith('#'):                     #---- skip the commented line
            continue
        ent         = ent.replace(' ', '')          #--- remove all white space
        atemp       = re.split(':', ent)
        key_val     = atemp[0].strip()
        first_list  = re.split('\,', atemp[1])      #--- database name
        second_list = re.split('\,', atemp[2])      #--- ocat data name
        test        = atemp[2].replace(' ', '')

        alist = [first_list]
#
#--- if the second list is empty, copy the first list in the position
#
        if test  == '':
            alist.append(first_list)
        else:
            alist.append(second_list)
        table_dict[key_val] = alist

    return table_dict

#----------------------------------------------------------------------------------------------------------------
#-- extract: extract value from the sql database                                                              ---
#----------------------------------------------------------------------------------------------------------------

def extract(database, table_dict, criteria, table_name):

    """
    extract value from the sql database
    Input:  database     database name
            table_dict   a dictionary containing a list of parameters
            criteria     data extraction criteria
            table_name   a name of the table list
    Output  temp_dict    a dictionary containing parameter <---> [value, value]
    """
#
#--- it is possible that the database name and Ocat Data page name are different
#
    name_list = table_dict[table_name][0]
    outnames  = table_dict[table_name][1]
#
#--- preparing to call SQL database
#
    cmd = 'select '
    for i in range(0, len(name_list)):
        cmd = cmd + name_list[i]
        if i < len(name_list) -1:
            cmd = cmd + ', '
    cmd = cmd + ' ' + criteria
#
#--- call SOL and get the result in a dictionary form
#
    out = readSQL(cmd, database)
#
#--- putting the result in the dictionary form
#
    t_dict = {}
    t_name = []
    for ent in outnames:
        name = ent.strip()              #--- clean up the key values: remove leading and trailing white space
        t_dict[name] = ['','']
        t_name.append(name)

    if out is not  None:
#
#---- value parts are in a list form. the first one will be "original" value and the second one is "current" value
#
        for i in range(0,len(name_list)):
            ent =  name_list[i]
            t_dict[t_name[i]] = [out[ent], out[ent]]

    return t_dict

#----------------------------------------------------------------------------------------------------------------
#-- extract_order_entry: extract value from the sql database: multiple order case                             ---
#----------------------------------------------------------------------------------------------------------------

def extract_order_entry(database, table_dict, criteria, table_name, oname):

    """
    extract value from the sql database: multiple order case
    Input:  database     database name
            table_dict   a dictionary containing a list of parameters
            criteria     data extraction criteria
            table_name   a name of the table list
    Output  temp_dict    a dictionary containing parameter <---> [[value for order=1, value for order=2...], [...]]

    """
    name_list = table_dict[table_name][0]
    outnames  = table_dict[table_name][1]

    for ent in outnames:
        exec "alist_%s = []" % (ent)
#
#--- loop around higher order: if fails at ordr= <ordr>, that will be the ordr of this value
#
    ordr = 0
    for i in range(1, 30):
        ocriteria = criteria + str(i)
        adict = extract(database, table_dict, ocriteria, table_name)

        chk = 0
        for ent in adict.values():
            if (ent[0] != '') and (ent[0] is not None):
                chk += 1
        if chk > 0:
            ordr = i
            for ent in outnames:
                exec "alist_%s.append(adict['%s'][0])" % (ent, ent)

    temp_dict = {oname : [ordr, ordr]}
    for ent in outnames:
        exec "toutnames = alist_%s" % (ent)
        tlist2 = copy.deepcopy(toutnames)
        exec "temp_dict['%s'] = [toutnames, tlist2]" %(ent)

    return temp_dict

#----------------------------------------------------------------------------------------------------------------
#-- readSQL: or a given sql command and database, return result in dictionary form                            ---
#----------------------------------------------------------------------------------------------------------------

def readSQL(cmd, database, all=''):

    """
    for a given sql command and database, return result in dictionary form.

	sqlbeta, sybase 15.5, on ascarc7/solaris
	arcbeta, sybase 15.7, on larcbeta/RHEL6.5


    """

    db_user = 'mtaops_internal_web'
    server  = 'ocatsqlsrv'
#    server  = 'sqlbeta'

    line      = pass_dir + '.targpass_internal'
    f         = open(line, 'r')
    db_passwd = f.readline().strip()
    f.close()

    db = dbi.DBI(dbi='sybase', server=server, user=db_user, passwd=db_passwd, database=database)

    if all == '':
        outdata = db.fetchone(cmd)
    else:
        outdata = db.fetchall(cmd)

    return outdata


#----------------------------------------------------------------------------------------------------------------
#--- find_monitor_obs: find obsid on a monitor list                                                           ---
#----------------------------------------------------------------------------------------------------------------

def find_monitor_obs(obsid):

    """
    for a given obsid, check all other obsid on the same monitor list. return monitor_list.

    """
    cmd_p = 'select obsid from target where pre_id='
    cmd_o = 'select pre_id from target where obsid='

    m_list = []
    mlist = find_series(obsid, m_list, cmd_p, cmd_o)
    mlist = find_series(obsid, m_list, cmd_o, cmd_p)


    m_list.sort()
    temp_list = set(m_list)
    monitor_list = list(temp_list)

    return monitor_list
    
#----------------------------------------------------------------------------------------------------------------
#-- find_series: for given obsid, database, and monitor list, exract series of pre-id realated to obsid       ---
#----------------------------------------------------------------------------------------------------------------

def find_series(obsid, m_list, cmd1, cmd2):

    """
    for given obsid, database, and monitor list, exract series of pre-id realated to obsid. return monitor_list
    this one checks obsid in decreasing manner.

    """

    m_list.append(obsid)
    p_list = []

    try:
        cmd =  cmd1 + str(obsid)
        out = readSQL(cmd, 'axafocat', all='all')

        ecnt = len(out)
        if ecnt == 0:
            return m_list
    except:
        return m_list

    for i in range(0, ecnt):
        obs = out[i][0]           #--- obs is obsid
        if chkNumeric(obs) == False:
            continue

        p_list.append(obs)

        try:
            cmd =  cmd2 + str(obs)
            out2 = readSQL(cmd, 'axafocat', all ='all')
            ecnt2 = len(out2)
            for j in range(0, ecnt2):
                obs2 = out2[j][0]
                if chkNumeric(obs2) == False:
                    continue

                if obs2 != obsid:
                    m_list = find_series(obs2, m_list, cmd2, cmd1)
        except:
            pass

    skip = 0
    for ent in m_list:
        for comp in  p_list:
            if ent == comp:
                skip = 1
                break
        if skip == 1:
            break

    if skip == 0:
        for ent in p_list:
            m_list = find_series(ent, m_list, cmd1, cmd2)

    return m_list

#-----------------------------------------------------------------------------------------------------------------
#--- chkNumeric: checkin entry is numeric value                                                                ---
#-----------------------------------------------------------------------------------------------------------------

def chkNumeric(elm):

    """
    check the entry is numeric. If so return True, else False.
    """

    try:
        test = float(elm)
    except:
        return False
    else:
        return True


#----------------------------------------------------------------------------------------------------------------
#-- find_planned_roll: get planned roll from mp web page                                                      ---
#----------------------------------------------------------------------------------------------------------------

def find_planned_roll(obsid):

    """
    get planned roll from mp web page
    Input:  obsid        --- obsid
            it also uses information from <obss_ss>/mp_long_term file
    Output: planned_roll --- mp planned roll angle
    """
    file = obs_ss + '/mp_long_term'
    f    = open(file, 'r')
    data = [line.strip() for line in f.readlines()]
    f.close()

    planned_roll = ''
    for ent in data:
        atemp = re.split(':', ent)

        if atemp[0] == obsid:
            planned_roll = atemp[1]
            break

    return planned_roll

#----------------------------------------------------------------------------------------------------------------
#-- TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST  ---
#----------------------------------------------------------------------------------------------------------------

class TestFunctions(unittest.TestCase):
    """
    testing functions
    """

#------------------------------------------------------------

    def test_param_table(self):
        table_dict = read_db_param_table()
        self.assertEquals(table_dict['flag_table'][0],['phase_constraint_flag', 'dither_flag', 'roll_flag', 'window_flag', 'spwindow_flag'])

#------------------------------------------------------------

    def test_extract(self):
        table_dict = read_db_param_table()
        obsid      = 14141
        database   = 'axafocat'
        criteria   = ' from target where obsid=' + str(obsid)
        table_name = 'target_table'
        odict      = extract(database, table_dict, criteria, table_name)
        self.assertEquals(odict['z_det_offset'], [-0.25, -0.25])
        self.assertEquals(odict['tooid'], [936, 936])
        self.assertEquals(odict['type'], ['TOO', 'TOO'])
        self.assertEquals(odict['ra'], [162.53708333333333, 162.53708333333333])
        self.assertEquals(odict['forder_cnt_rate'], [None, None])

#------------------------------------------------------------

    def test_monitor(self):
        alist     = [14330, 12612, 12613, 12614, 12615]
        test_list = find_monitor_obs(12614)
        self.assertEquals(alist, test_list)
 
#        test = find_monitor_obs(15442)
#        print "I AM HERE: " + str(test)

#------------------------------------------------------------

    def test_order(self):
        obsid      = 12109
        database   = 'axafocat'
        table_dict = read_db_param_table()
        criteria   = ' from timereq where  obsid=' + str(obsid) + ' and ordr='
        table_name = 'time_table'
        oname      = 'time_ordr'

        out        = extract_order_entry(database, table_dict, criteria, table_name, oname)
        self.assertEquals(str(out['tstop'][1][0]),'Feb 28 2010 11:59PM')

#------------------------------------------------------------

class TestClassMethods(unittest.TestCase):
    """
    testing class methods
    """

    def test_Values(self):
        db = OcatDB(12109)

        list = db.paramList() 
        self.assertEquals(list[1], 'window_constraint')

        val = db.origValue('si_mode')
        self.assertEquals(val, 'CC_0010C')
        val = db.newValue('si_mode')
        self.assertEquals(val, 'CC_0010C')

        db.updateValue('si_mode', 'xxxxxx')
        val = db.newValue('si_mode')
        self.assertEquals(val, 'xxxxxx')

        db.addParam('newparam', '1234')
        val = db.newValue('newparam')
        self.assertEquals(val, '1234')

        db.addNewOrder('roll_ordr', ['P', 'Y', '90', '14'])
        val = db.newValue('roll')
        self.assertEquals(val[0], '90')

        test = db.getTableEntries('dither_table')
        self.assertEquals(test[0], 'y_amp')


#---------------------------------------------------------------------------------------------

if __name__ == '__main__':

    """
    if you just run this script, it will run a test mode.
    """
    unittest.main()
