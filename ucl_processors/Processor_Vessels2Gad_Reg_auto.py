"""Processor associated to Spider_Vessels2Gad_Reg.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Processor name: Processor_Vessels2Gad_Reg
Creation date:  2017-05-03 10:52:08.743152
Purpose:        Run the Vessesl to Gad registration.
"""

# Python packages import
import os
import logging
from dax import XnatUtils, SessionProcessor

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Run the Vessesl to Gad registration."
__processor_name__ = "Processor_Vessels2Gad_Reg"
__modifications__ = """2017-05-03 10:52:08.743152 - Original write"""

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
HOME = os.path.expanduser("~")
DEFAULT_SPIDER_PATH = os.path.join(HOME, 'Xnat-management/ucl_processing/\
ucl_spiders/Spider_Vessels2Gad_Reg_v1_0_0.py')
DEFAULT_WALLTIME = '02:00:00'
DEFAULT_MEM = 4048
DEFAULT_PPN = 4
DEFAULT_WORKING_DIR = '/scratch0/dax/'
SCAN_RESOURCE = 'NIFTI'
GAD_SCAN = ['1mm Ax FSPG + Gad']
VESSELS_SCANS = ['Ax Inhance 3D MR*']
NIPYPE_EXE = 'perform_vessels2gad_registration.py'

# Format for the spider command line
SPIDER_FORMAT = """python {spider} \
-a {assr} \
-d {dir} \
--exe {exe} \
--gad {gad} \
--vessels {vessels} \
--vesselsids {vesselsids} \
--suffix "{suffix_proc}" \
--openmp_core {number_core} \
--working_dir "{working_dir}"
"""


class Processor_Vessels2Gad_Reg(SessionProcessor):
    """Processor class for Vessels2Gad_Reg that runs on a session.

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider
    #
    # ADD MORE PARAMETERS AS NEEDED HERE AND IN THE __INIT__
    #
    :param suffix: suffix to the spider
    """

    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 gad_type=GAD_SCAN, vessels_type=VESSELS_SCANS,
                 exe=NIPYPE_EXE, walltime=DEFAULT_WALLTIME,
                 mem_mb=DEFAULT_MEM, ppn=DEFAULT_PPN,
                 suffix_proc='', working_dir=DEFAULT_WORKING_DIR):
        """Entry point for Processor_Vessels2Gad_Reg Class."""
        super(Processor_Vessels2Gad_Reg,
              self).__init__(walltime, mem_mb, spider_path, version,
                             suffix_proc=suffix_proc)
        self.gad_type = gad_type
        self.vessels_type = vessels_type
        self.exe = exe
        self.ppn = ppn
        self.working_dir = working_dir

    def has_inputs(self, csess):
        """Method overridden from base class.

        By definition:
            status = 0  -> NEED_INPUTS,
            status = 1  -> NEED_TO_RUN
            status = -1 -> NO_DATA
            qcstatus needs a value only when -1 or 0.
        You need to set qcstatus to a short string that explain
        why it's no ready to run. e.g: No NIFTI

        :param csess: object csess define in dax.XnatUtils
                      (see XnatUtils in dax for information)
        :return: status, qcstatus
        """
        gad_cscans = XnatUtils.get_good_cscans(csess, self.gad_type)
        if not gad_cscans:
            LOGGER.debug('Vessel Registration: No GAD scan found.')
            return -1, 'No GAD found'
        elif len(gad_cscans) > 1:
            LOGGER.debug('Vessel Registration: Too many GAD scans found.')
            return 0, 'Too many GAD found'

        vessels_cscans = XnatUtils.get_good_cscans(csess, self.vessels_type)
        if not vessels_cscans:
            LOGGER.debug('Vessel Registration: No Vessels scans found.')
            return -1, 'No Vessels found'

        return 1, None

    def get_file_path(self, cscans, resource):
        """Method to get the file path on XNAT for the scans

        :param csess: object csess define in dax.XnatUtils
                      (see XnatUtils in dax for information)
        :return: list of paths
        """
        filepaths = list()
        path_tmp = 'xnat:/project/{0}/subject/{1}/experiment/{2}/scan/{3}/\
resource/{4}'

        for cscan in cscans:
            scan_info = cscan.info()
            filepaths.append(path_tmp.format(scan_info['project_id'],
                                             scan_info['subject_label'],
                                             scan_info['session_label'],
                                             scan_info['id'],
                                             resource))
        return filepaths

    def get_cmds(self, assessor, jobdir):
        """Method to generate the spider command for cluster job.

        :param assessor: pyxnat assessor object
        :param jobdir: jobdir where the job's output will be generated
        :return: command to execute the spider in the job script
        """
        assessor_label = assessor.label()
        proj_label = assessor.parent().parent().parent().label()
        subj_label = assessor.parent().parent().label()
        sess_label = assessor.parent().label()

        csess = XnatUtils.CachedImageSession(assessor._intf, proj_label,
                                             subj_label, sess_label)

        gad_cscans = XnatUtils.get_good_cscans(csess, self.gad_type)
        gad = self.get_file_path(gad_cscans[0], SCAN_RESOURCE)

        vessels_cscans = XnatUtils.get_good_cscans(csess, self.vessels_type)
        vessels = self.get_file_path(vessels_cscans, SCAN_RESOURCE)
        vessel_ids = [ves.info()['ID'] for ves in vessels_cscans]

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   assr=assessor_label,
                                   dir=jobdir,
                                   exe=self.exe,
                                   gad=','.join(gad),
                                   vessels=','.join(vessels),
                                   vesselsids=','.join(vessel_ids),
                                   number_core=self.ppn,
                                   working_dir=self.working_dir)

        return [cmd]
