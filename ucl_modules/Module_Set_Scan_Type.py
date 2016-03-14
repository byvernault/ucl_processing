""" Module to set the scan type from a text file """
from dax import SessionModule
import os
import logging
import fileinput

LOGGER = logging.getLogger('dax')

DEFAULT_TPM_PATH = 'NONE'
DEFAULT_MODULE_NAME = 'Set_Scan_Type'
DEFAULT_TEXT_REPORT = 'ERROR/WARNING for setting scan type:\n'
DEFAULT_SCANTYPE_FILE_PATH = ''
RESOURCE_FLAG_NAME = 'Module_Set_Scan_Type'

class Module_Set_Scan_Type(SessionModule):
    """ Module to set the scan types on XNAT """
    def __init__(self, mod_name=DEFAULT_MODULE_NAME, directory=DEFAULT_TPM_PATH, email=None, text_report=DEFAULT_TEXT_REPORT,
                 scantype_file=DEFAULT_SCANTYPE_FILE_PATH):
        """ init function overridden from base-class"""
        super(Module_Set_Scan_Type, self).__init__(mod_name, directory, email, text_report=text_report)
        self.scantype_file = scantype_file
        self.old_scan_type_list = list()
        self.new_scan_type_list = list()

    def prerun(self, settings_filename=''):
        """ prerun function overridden from base-class"""
        #read text file
        if os.path.exists(self.scantype_file):
            for line in fileinput.input(self.scantype_file):
                stringline = line.strip().split(';')
                scans_type = stringline[0].strip().split('=')
                if len(scans_type) > 1:
                    self.old_scan_type_list.append(scans_type[0])
                    self.new_scan_type_list.append(scans_type[1])

    def afterrun(self, xnat, project):
        """ afterrun function overridden from base-class"""
        #send report
        if self.send_an_email:
            self.send_report()

    def needs_run(self, csess, xnat):
        """ needs_run function overridden from base-class
                csession = CacheSession object from XnatUtils
            return True or False
        """
        return self.has_flag_resource(csess, RESOURCE_FLAG_NAME)

    def run(self, session_info, session_obj):
        """ main function to run on the session """
        LOGGER.debug('Setting Scan type')
        for scan_obj in session_obj.scans().fetchall('obj'):
            if scan_obj.attrs.get('xnat:imageScanData/type') in self.old_scan_type_list:
                #INDEX
                index_scantype = self.old_scan_type_list.index(scan_obj.attrs.get('xnat:imageScanData/type'))
                #set the new scan type
                scan_obj.attrs.set('xnat:imageScanData/type', self.new_scan_type_list[index_scantype])

        session_obj.resource(RESOURCE_FLAG_NAME).create()
