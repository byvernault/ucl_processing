"""Processor associated to Spider_ADC_MAP.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Processor name: Processor_ADC_MAP
Creation date:  2016-08-18 13:14:30.689905
Purpose:        Generate ADC map from Verdict registered scans
"""

# Python packages import
import os
import logging
from dax import XnatUtils, ScanProcessor

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Generate ADC map from Verdict registered scans"
__processor_name__ = "Processor_Verdict"
__modifications__ = """2016-08-18 13:14:30.689905 - Original write"""

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
HOME = os.path.expanduser("~")
DEFAULT_SPIDER_PATH = os.path.join(HOME, 'Xnat-management/ucl_processing/\
ucl_spiders/', 'Spider_ADC_MAP_v1_0_0.py')
DEFAULT_WALLTIME = '00:30:00'
DEFAULT_MEM = 6048
DEFAULT_MATLAB_CODE = os.path.join(HOME, 'Code', 'matlab')
DEFAULT_SCAN_TYPES = ['ep2d_DWI_ax_high_res$']

# Format for the spider command line
SPIDER_FORMAT = """python {spider} \
-p {proj} \
-s {subj} \
-e {sess} \
-c {scan} \
-d {dir} \
--suffix "{suffix_proc}" \
--mc {mc} \
"""


class Processor_ADC_MAP(ScanProcessor):
    """Processor class for Verdict that runs on a session.

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider

    :param suffix: suffix to the spider
    """

    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 scan_types=DEFAULT_SCAN_TYPES,
                 matlab_code=DEFAULT_MATLAB_CODE,
                 walltime=DEFAULT_WALLTIME, mem_mb=DEFAULT_MEM,
                 suffix_proc=''):
        """Entry point for Processor_Verdict Class."""
        super(Processor_ADC_MAP,
              self).__init__(scan_types, walltime, mem_mb, spider_path,
                             version, suffix_proc=suffix_proc)
        self.mc = matlab_code

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
        if XnatUtils.has_resource(cscan, 'DICOM'):
            return 1, None

        LOGGER.debug('ADC MAP: DICOM not found.')
        return 0, 'No DICOM'

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
                                   dir=jobdir,
                                   mc=self.mc,
                                   suffix_proc=self.suffix_proc)

        return [cmd]
