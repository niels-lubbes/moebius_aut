'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 17, 2017
@author: Niels Lubbes
'''

from moebius_aut.class_ma_tools import MATools

class TestClassMATools:


    def test__p( self ):

        ma = MATools()

        ma.filter( None )
        assert ma.p( 'Hello world!' ) != None

        ma.filter( 'another_class.py' )
        assert ma.p( 'No output since called from another class.' ) == None

        ma.filter_unset()
        assert ma.p( 'Filter is disabled so output this string.' ) != None

        ma.filter_reset()
        assert ma.p( 'Filter is enabled again so do not output.' ) == None

        ma.filter( 'test_class_ma_tools.py' )
        assert ma.p( 'Only output if called from this class' ) != None


    def test__tool_dct( self ):

        ma = MATools()
        ma2 = MATools()

        # watch out to not use the default file name
        # otherwise it might take long to load the data
        test_fname = 'test_tools'
        key = 'test__tool_dct'

        dct = ma.get_tool_dct( fname = test_fname )
        dct[key] = True
        ma.save_tool_dct( fname = test_fname )

        assert key in ma.get_tool_dct( fname = test_fname )
        assert key in ma2.get_tool_dct( fname = test_fname )

        ma.set_enable_tool_dct( False )
        assert key not in ma.get_tool_dct( fname = test_fname )
        assert key not in ma2.get_tool_dct( fname = test_fname )

        ma.set_enable_tool_dct( True )
        assert key in ma.get_tool_dct( fname = test_fname )
        assert key in ma2.get_tool_dct( fname = test_fname )


