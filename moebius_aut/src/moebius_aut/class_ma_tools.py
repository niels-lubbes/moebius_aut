'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 17, 2017
@author: Niels Lubbes
'''
from sage.all import *

import inspect
import time
import sys


class MATools():
    '''
    For accessing static variables in python see for example:
    <http://stackoverflow.com/questions/68645/static-class-variables-in-python>    
    '''

    # Private dictionary object for caching result
    # used by ".get_tool_dct()" and ".save_tool_dct()".
    # If "enable_tool_dct" is false then caching in
    # disabled. This is useful for example in test
    # methods. However, it should be noted that it
    # could take a long time to compute the data.
    #
    __tool_dct = None
    __enable_tool_dct = True

    # private variable for timer
    #
    __start_time = None
    __end_time = None

    # private static variables used by ".verbose()"
    #
    __filter_fname = None
    __prev_filter_fname = None


    @staticmethod
    def filter( filter_fname ):
        '''
        It is adviced to access this method 
        as MATools.filter().  
        
        INPUT:
            - "filter_fname" -- File name 
        '''
        MATools.__filter_fname = filter_fname
        MATools.__prev_filter_fname = filter_fname

    @staticmethod
    def filter_unset():
        '''
        Output via ".out" will not be surpressed.
        '''
        MATools.__filter_fname = None

    @staticmethod
    def filter_reset():
        '''
        Resets filter state to before previous ".filter_unset()" call.
        '''
        MATools.__filter_fname = MATools.__prev_filter_fname

    @staticmethod
    def p( *arg_lst ):
        '''
        INPUT:
            - "*arg_lst" -- List of arguments.
        OUTPUT:
            - If ".filter_on(<fname>)" has been called and the file name
              of the calling module does not coincide with <fname>
              then the output is surpressed and "None" is returned.
                            
              Otherwise, this method prints arguments to "sys.stdout" 
              together with reflection info from "inspect.stack()".
              Additional returns the output string.
              
              Call ".filter_off()" to turn off filter, such that
              all output is send to "sys.stdout".  
                                   
        '''
        # collect relevant info from stack trace
        sk_lst_lst = inspect.stack()
        file_name = str( sk_lst_lst[1][1] )
        line = str( sk_lst_lst[1][2] )
        method_name = str( sk_lst_lst[1][3] )

        # only output when op is called from "op.input_file_name"
        if MATools.__filter_fname != None:
            if not file_name.endswith( MATools.__filter_fname ):
                return

        # construct output string
        s = method_name + '(' + line + ')' + ': '
        for arg in arg_lst:
            s += str( arg ) + ' '

        # print output
        print( s )
        sys.stdout.flush()

        return s


    @staticmethod
    def set_enable_tool_dct( enable_tool_dct ):
        MATools.__enable_tool_dct = enable_tool_dct

    @staticmethod
    def get_tool_dct( fname = 'ma_tools' ):
        '''
        INPUT:
            - "fname" -- Name of file without extension.
        OUTPUT:
            - Sets static private variable "__tool_dct" 
              in memory from file "<local path>/<fname>.sobj"
              if called for the first time.
              
            - Returns ".__tool_dct" if ".__enable_tool_dct==True" 
              and "{}" otherwise.
        '''
        if not MATools.__enable_tool_dct:
            MATools.filter_unset()
            MATools.p( 'Caching is disabled!' )
            MATools.filter_reset()
            return {}

        path = os.path.dirname( os.path.abspath( __file__ ) ) + '/'
        file_name = path + fname
        if MATools.__tool_dct == None:

            MATools.filter_unset()
            try:

                MATools.p( 'Loading from:', file_name )
                MATools.__tool_dct = load( file_name )

            except Exception as e:

                MATools.p( 'Cannot load ".__tool_dct": ', e )
                MATools.__tool_dct = {}

            MATools.filter_reset()

        return MATools.__tool_dct


    @staticmethod
    def save_tool_dct( fname = 'ma_tools' ):
        '''
        INPUT:
            - "fname" -- Name of file without extension.        
        OUTPUT:
            - Saves ".__tool_dct" to  "fname" if ".enable_tool_dct==True" 
              otherwise do nothing.
        '''
        if not MATools.__enable_tool_dct:
            MATools.filter_unset()
            MATools.p( 'Data is not saved to disk!' )
            MATools.filter_reset()
            return

        path = os.path.dirname( os.path.abspath( __file__ ) ) + '/'
        file_name = path + fname

        MATools.filter_unset()
        MATools.p( 'Saving to:', file_name )
        MATools.filter_reset()

        save( MATools.__tool_dct, file_name )


    @staticmethod
    def start_timer():
        '''
        OUTPUT:
            - Prints the current time and starts timer.
        '''
        # get time
        MATools.__start_time = time.clock()  # set static variable.

        MATools.filter_unset()
        MATools.p( 'start time =', MATools.__start_time )
        MATools.filter_reset()


    @staticmethod
    def end_timer():
        '''
        OUTPUT:
            - Prints time passed since last call of ".start_timer()".
        '''
        MATools.__end_time = time.clock()
        passed_time = MATools.__end_time - MATools.__start_time

        MATools.filter_unset()
        MATools.p( 'time passed =', passed_time )
        MATools.filter_reset()




