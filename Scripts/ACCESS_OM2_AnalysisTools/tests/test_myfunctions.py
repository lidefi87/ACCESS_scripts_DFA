from ACCESS_Analysis import myfunctions

def test_getACCESSdata():
    assert myfunctions.getACCESSdata('aice', '1958-02', '1959-03', freq = '1 daily', 
                                     ses = cc.database.create_session(), exp = '01deg_jra55v140_iaf_cycle2',
                                    ice_data = True)

# SO = zsf.getACCESSdata(varInt, stime[i], etime[i], freq = freq, ses = session,\
#                            minlat = -90, maxlat = -45, exp = exp, ice_data = True)