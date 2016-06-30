"""Processor associated to Spider_GIF_Parcellation.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Processor name: Processor_GIF_Parcellation
Creation date:  2016-02-08 15:17:05.774006
Purpose:        Processor for Parcellation of the brain using GIF:
                Geodesic Information Flow
"""

# Python packages import
import os
import logging
from dax import XnatUtils, ScanProcessor

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Processor for Parcellation of the brain using GIF: \
Geodesic Information Flow"
__processor_name__ = "Processor_GIF_Parcellation"
__modifications__ = "2016-03-15 14:56 - Adding working_dir options"

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
HOME = os.path.expanduser("~")
DEFAULT_SPIDER_PATH = os.path.join(HOME, 'Xnat-management/ucl_processing/\
ucl_spiders/Spider_GIF_Parcellation_v1_0_0.py')
DEFAULT_WALLTIME = '48:00:00'
DEFAULT_MEM = 3850
DEFAULT_PPN = 4
DEFAULT_WORKING_DIR = '/scratch0/dax/'
DEFAULT_TEMPLATE = '/cluster/project0/GIF/template-database-r2.1/db.xml'
DEFAULT_GIF_PATH = 'perform_gif_propagation.py'
DEFAULT_SCAN_TYPES = ['T1', 'MPRAGE']  # ADD SCAN TYPES

# Format for the spider command line
SPIDER_FORMAT = """python {spider} \
-p {proj} \
-s {subj} \
-e {sess} \
-c {scan} \
-d {dir} \
--suffix "{suffix_proc}" \
--dbt {template} \
--gif {gif_path} \
--openmp_core {number_core} \
--working_dir "{working_dir}"
"""


class Processor_GIF_Parcellation(ScanProcessor):
    """Processor class for GIF_Parcellation that runs on a scan.

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider
    :param ppn: number of core used
    :param scan_types: scan types on XNAT that the spider should run on
    :param gif: executable path for perform_gif_propagation
    :patam db_template: db.xml file containing the template for atlases
    :param suffix: suffix to the spider
    :param working_dir: working directory for spider temp files
    """

    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 walltime=DEFAULT_WALLTIME, mem_mb=DEFAULT_MEM,
                 ppn=DEFAULT_PPN, db_template=DEFAULT_TEMPLATE,
                 gif=DEFAULT_GIF_PATH, scan_types=DEFAULT_SCAN_TYPES,
                 suffix_proc='',  working_dir=DEFAULT_WORKING_DIR):
        """Entry point for Processor_GIF_Parcellation Class."""
        super(Processor_GIF_Parcellation,
              self).__init__(scan_types, walltime, mem_mb, spider_path,
                             version, ppn=ppn, suffix_proc=suffix_proc)
        self.db_template = db_template
        self.gif = gif
        self.working_dir = working_dir

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

        LOGGER.debug('GIF Parcellation: NIFTI not found.')
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

        working_dir = os.path.join(self.working_dir, assr_label)

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   proj=proj_label,
                                   subj=subj_label,
                                   sess=sess_label,
                                   scan=scan_label,
                                   dir=jobdir,
                                   suffix_proc=self.suffix_proc,
                                   template=self.db_template,
                                   gif_path=self.gif,
                                   number_core=self.ppn,
                                   working_dir=working_dir)

        return [cmd]
