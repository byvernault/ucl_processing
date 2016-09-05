"""Processor associated to Spider_Sample_GM_Segment.

Author:         Benjamin Yvernault & Dave Cash
contact:        b.yvernault@ucl.ac.uk
Processor name: Processor_Sample_GM_Segment
Creation date:  2016-09-05 10:25:39.292368
Purpose:        Run Sample GM Segment spider
"""

# Python packages import
import os
import logging
from dax import XnatUtils, ScanProcessor

__author__ = "Benjamin Yvernault & Dave Cash"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Run Sample GM Segment spider"
__processor_name__ = "Processor_Sample_GM_Segment"
__modifications__ = """2016-09-05 10:25:39.292368 - Original write"""

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
HOME = os.path.expanduser("~")
DEFAULT_SPIDER_PATH = os.path.join(HOME, 'Xnat-management/ucl_processing/\
ucl_spiders/', 'Spider_Sample_GM_Segment_v1_0_0.py')
DEFAULT_WALLTIME = '01:00:00'
DEFAULT_MEM = 3048

DEFAULT_SCAN_TYPES = []  # ADD SCAN TYPES
DEFAULT_SPM12 = os.path.join(HOME, 'Code', 'spm12')
DEFAULT_MATLAB_CODE = os.path.join(HOME, 'Code', 'matlab')

# Format for the spider command line
SPIDER_FORMAT = """python {spider} \
-p {proj} \
-s {subj} \
-e {sess} \
-c {scan} \
-d {dir} \
--spm12 {spm12} \
--mc {matlab_code} \
--suffix "{suffix_proc}"
"""


class Processor_Sample_GM_Segment(ScanProcessor):
    """Processor class for Sample_GM_Segment that runs on a scan.

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider
    :param scan_types: scan types on XNAT that the spider should run on
    :param matlab_code: path to matlab function Sample_GM_Segment
    :patam spm12: path to spm12
    :param suffix: suffix to the spider
    """

    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 walltime=DEFAULT_WALLTIME, mem_mb=DEFAULT_MEM,
                 matlab_code=DEFAULT_MATLAB_CODE, spm12=DEFAULT_SPM12,
                 scan_types=DEFAULT_SCAN_TYPES, suffix_proc=''):
        """Entry point for Processor_Sample_GM_Segment Class."""
        super(Processor_Sample_GM_Segment,
              self).__init__(scan_types, walltime, mem_mb, spider_path,
                             version, suffix_proc=suffix_proc)
        self.matlab_code = matlab_code
        self.spm12 = spm12

    def has_inputs(self, cscan):
        """Method overridden from base class.

        By definition:
            status = 0  -> NEED_INPUTS,
            status = 1  -> NEED_TO_RUN
            status = -1 -> NO_DATA
            qcstatus needs a value only when -1 or 0.
        You need to set qcstatus to a short string that explain
        why it's no ready to run. e.g: No NIFTI

        :param cscan: object cscan define in dax.XnatUtils
                      (see XnatUtils in dax for information)
        :return: status, qcstatus
        """
        if XnatUtils.is_cscan_unusable(cscan):
            return -1, 'Scan unusable'

        # Check has_resource PARREC or NIFTI
        if XnatUtils.has_resource(cscan, 'NIFTI'):
            return 1, None

        LOGGER.debug('Sample GM Segemnt: NIFTI not found.')
        return 0, 'No NIFTI'

    def get_cmds(self, assessor, jobdir):
        """Method to generate the spider command for cluster job.

        :param assessor: pyxnat assessor object
        :param jobdir: jobdir where the job's output will be generated
        :return: command to execute the spider in the job script
        """
        proj_label = assessor.parent().parent().parent().label()
        subj_label = assessor.parent().parent().label()
        sess_label = assessor.parent().label()
        assr_label = assessor.label()
        scan_label = assr_label.split('-x-')[3]

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   proj=proj_label,
                                   subj=subj_label,
                                   sess=sess_label,
                                   scan=scan_label,
                                   matlab_code=self.matlab_code,
                                   spm12=self.spm12,
                                   dir=jobdir,
                                   suffix_proc=self.suffix_proc)

        return [cmd]
