'''
    Author:         Benjamin Yvernault
    contact:        b.yvernault@ucl.ac.uk
    Processor name: Processor_Registration2Ref
    Creation date:  2016-02-24 14:06:18.168057
    Purpose:        Register scans from a session with a specific types to one reference scan (same session)
'''

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Register scans from a session with a specific types to one reference scan (same session)"
__processor_name__ = "Processor_Registration2Ref"
__modifications__ = "2016-02-24 14:06:18.168057 - Original write"

# Python packages import
import os
import logging
from dax import XnatUtils, SessionProcessor

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
HOME = os.path.expanduser("~")
DEFAULT_SPIDER_PATH = os.path.join(HOME, 'Xnat-management/ucl_processing/ucl_spiders/Spider_Registration2Ref_v1_0_0.py')
DEFAULT_WALLTIME = '00:30:00'
DEFAULT_MEM = 2048
DEFAULT_PPN = 1
DEFAULT_TYPES = 'all'
DEFAULT_REG_ALADIN = '/share/apps/cmic/niftypipe_deps/bin/reg_aladin'

# Format for the spider command line
SPIDER_FORMAT = '''python {spider} -p {proj} -s {subj} -e {sess} -d {dir} --suffix "{suffix_proc}" --scansID {sources} --scanRef {target} --regAladin {regaladin} --openm_core {number_core}'''

class Registration2Ref_Processor(SessionProcessor):
    '''
    Processor class for Registration2Ref that runs on a session

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider
    :param target_type: scan type for the reference
    :param sources_type: scan types to register to the reference
    :param reg_aladin_exe: path to reg_aladin executable
    :param suffix: suffix to the spider
    '''
    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 walltime=DEFAULT_WALLTIME, mem_mb=DEFAULT_MEM, ppn=DEFAULT_PPN,
                 target_type=None, sources_type=DEFAULT_TYPES,
                 reg_aladin_exe=DEFAULT_REG_ALADIN, suffix_proc=''):
        super(Registration2Ref_Processor, self).__init__(walltime, mem_mb, spider_path, version,
                                                         ppn=ppn, suffix_proc=suffix_proc)
        # Reference
        if not target_type:
            raise Exception("Processor_Registration2Ref: reference type not set. Please edit your settings file.")
        self.target_type = target_type
        # Images to register
        if isinstance(sources_type, list):
            self.sources_type = sources_type
        elif isinstance(sources_type, str):
            if sources_type != 'all':
                sources_type = sources_type.split(',')
        if not sources_type:
            raise Exception("Processor_Registration2Ref: scantypes to register not set. Please edit your settings file.")
        # Reg Aladin
        if not os.path.exists(reg_aladin_exe):
            raise Exception("Processor_Registration2Ref: reg_aladin path given not found: %s." % reg_aladin_exe)
        self.regaladin = reg_aladin_exe

    def has_inputs(self, csess):
        '''
        function overridden from base class
        By definition, status = 0 if still NEED_INPUTS, -1 if NO_DATA, 1 if NEED_TO_RUN
                       qcstatus needs a value only when -1 or 0.
        You can set qcstatus to a short string that explain why it's no ready to run.
            e.g: No NIFTI

        :param csess: object csess define in dax.XnatUtils (see XnatUtils in dax for information)
        :return: status, qcstatus
        '''
        # Check that there is only one scan usable with the reference type and that the reference file as a NIFTI
        target_cscans = XnatUtils.get_good_cscans(csess, self.target_type)
        if not target_cscans:
            LOGGER.debug('Processor_Registration2Ref: cannot run at all, no target image found')
            return -1, 'Target not found'
        if len(target_cscans) > 1:
            LOGGER.debug('Processor_Registration2Ref: cannot run at all, too many target images found')
            return 0, 'Too many target'
        if not XnatUtils.has_resource(target_cscans[0], 'NIFTI'):
            LOGGER.debug('Processor_Registration2Ref: cannot run, no NIFTI for target image')
            return 0, "no target's NIFTI"

        source_cscans = XnatUtils.get_good_cscans(csess, self.sources_type)
        if not source_cscans:
            LOGGER.debug('Processor_Registration2Ref: cannot run at all, no source image found')
            return -1, 'Sources not found'

        return 1, None

    def get_cmds(self, assessor, jobdir):
        '''
        function to generate the spider command for cluster job

        :param assessor: pyxnat assessor object
        :param jobdir: jobdir where the job's output will be generated
        :return: command to execute the spider in the job script
        '''
        proj_label = assessor.parent().parent().parent().label()
        subj_label = assessor.parent().parent().label()
        sess_label = assessor.parent().label()

        csess = XnatUtils.CachedImageSession(assessor._intf, proj_label, subj_label, sess_label)
        target_cscans = XnatUtils.get_good_cscans(csess, self.target_type)
        target_id = target_cscans[0].info()['ID']
        source_cscans = XnatUtils.get_good_cscans(csess, self.sources_type)
        sources_id = ','.join([sc.info()['ID'] for sc in source_cscans])

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   proj=proj_label,
                                   subj=subj_label,
                                   sess=sess_label,
                                   dir=jobdir,
                                   target=target_id,
                                   sources=sources_id,
                                   regaladin=self.regaladin,
                                   number_core=self.ppn,
                                   suffix_proc=self.suffix_proc)

        return [cmd]
