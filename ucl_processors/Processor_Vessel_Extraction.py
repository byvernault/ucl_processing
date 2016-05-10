"""Processor to extract vessel using niftkVesselExtractor.

Author:         Benjamin Yvernault
contact:        byvernault@gmail.com
Processor name: Processor_Vessel_Extraction
Creation date:  2016-05-10 11:28:41.384069
Purpose:        Extract the vessel from the T1/MPRAGE scan
"""

# Python packages import
import os
import logging
from dax import XnatUtils, ScanProcessor

__author__ = "Benjamin Yvernault"
__email__ = "byvernault@gmail.com"
__purpose__ = "Extract the vessel from the T1/MPRAGE scan"
__processor_name__ = "Processor_Vessel_Extraction"
__modifications__ = "2016-05-10 11:28:41.384069 - Original write"

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
HOME = os.path.expanduser("~")
CS_PATH = 'Xnat-management/ucl_processing/ucl_spiders'
DEFAULT_SPIDER_PATH = os.path.join(HOME, CS_PATH,
                                   'Spider_Vessel_Extraction_v1_0_0.py')
DEFAULT_WALLTIME = '02:00:00'
DEFAULT_MEM = 2048
DEFAULT_PIXEL_SIZE = '0.775438'
DEFAULT_SCAN_TYPES = ['T1', 'MPRAGE']  # ADD SCAN TYPES
DEFAULT_VESSEL_PATH = 'niftkVesselExtractor'

# Format for the spider command line
SPIDER_FORMAT = '''python {spider} -p {proj} -s {subj} -e {sess} -c {scan} \
--min {pixel_size} --vesselExtPath {vessel_path} -d {dir} \
--suffix "{suffix_proc}"'''


class Processor_Vessel_Extraction(ScanProcessor):
    """Processor class for Vessel_Extraction that runs on a scan.

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider
    :param scan_types: scan types on XNAT that the spider should run on
    :param vessel_path: path to niftkVesselExtractor
    :param suffix: suffix to the spider
    """

    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 walltime=DEFAULT_WALLTIME, mem_mb=DEFAULT_MEM,
                 scan_types=DEFAULT_SCAN_TYPES,
                 vessel_path=DEFAULT_VESSEL_PATH,
                 suffix_proc=''):
        """Init method."""
        super(Processor_Vessel_Extraction, self).__init__(scan_types,
                                                          walltime,
                                                          mem_mb,
                                                          spider_path,
                                                          version,
                                                          suffix_proc=suffix_proc)
        self.vessel_path = vessel_path

    def has_inputs(self, cscan):
        """function overridden from base class.

        By definition, status = 0 if still NEED_INPUTS, -1 if NO_DATA,
            1 if NEED_TO_RUN - qcstatus needs a value only when -1 or 0.
        You can set qcstatus to a short string that explain why it's not
        ready to run. e.g: No NIFTI

        :param cscan: object cscan define in dax.XnatUtils
         (see XnatUtils in dax for information)
        :return: status, qcstatus
        """
        if XnatUtils.is_cscan_unusable(cscan):
            return -1, 'Scan unusable'

        # Check has_resource PARREC or NIFTI
        if XnatUtils.has_resource(cscan, 'NIFTI'):
            return 1, None

        LOGGER.debug('Vessel Extraction: NIFTI not found.')
        return 0, 'No NIFTI'

    def get_cmds(self, assessor, jobdir):
        """function to generate the spider command for cluster job.

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
                                   vessel_path=self.vessel_path,
                                   dir=jobdir,
                                   suffix_proc=self.suffix_proc)

        return [cmd]
