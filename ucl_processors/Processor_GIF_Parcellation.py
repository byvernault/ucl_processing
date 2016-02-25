'''
    Author:         Benjamin Yvernault
    contact:        b.yvernault@ucl.ac.uk
    Processor name: Processor_GIF_Parcellation
    Creation date:  2016-02-08 15:17:05.774006
    Purpose:        Processor for Parcellation of the brain using GIF: Geodesic Information Flow
'''

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Processor for Parcellation of the brain using GIF: Geodesic Information Flow"
__processor_name__ = "Processor_GIF_Parcellation"
__modifications__ = "2016-02-08 15:17:05.774006 - Original write"

# Python packages import
import logging
from dax import XnatUtils, ScanProcessor

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
DEFAULT_SPIDER_PATH = '/home/byvernau/Xnat-management/spiders/Spider_GIF_Parcellation_v1_0_0.py'
DEFAULT_WALLTIME = '48:00:00'
DEFAULT_MEM = 3850
DEFAULT_PPN = 4
DEFAULT_TEMPLATE = '/cluster/project0/GIF/template-database-r2.1/db.xml'
DEFAULT_GIF_PATH = '/share/apps/cmic/niftypipe/bin/perform_gif_propagation.py'
DEFAULT_SCAN_TYPES = ['T1', 'MPRAGE'] # ADD SCAN TYPES

# Format for the spider command line
SPIDER_FORMAT = '''python {spider} -p {proj} -s {subj} -e {sess} -c {scan} -d {dir} --host "{host}" --suffix "{suffix_proc}" --dbt {template} --gif {gif_path}'''

class Processor_GIF_Parcellation(ScanProcessor):
    '''
    Processor class for GIF_Parcellation that runs on a scan

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider
    :param ppn: number of core used
    :param scan_types: scan types on XNAT that the spider should run on
    :param suffix: suffix to the spider
    '''
    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 walltime=DEFAULT_WALLTIME, mem_mb=DEFAULT_MEM, ppn=DEFAULT_PPN,
                 db_template=DEFAULT_TEMPLATE, gif=DEFAULT_GIF_PATH, xnat_host='',
                 scan_types=DEFAULT_SCAN_TYPES, suffix_proc=''):
        super(Processor_GIF_Parcellation, self).__init__(scan_types, walltime, mem_mb, spider_path,
                                                         version, ppn=ppn, suffix_proc=suffix_proc)
        self.xnat_host = xnat_host
        self.db_template = db_template
        self.gif = gif

    def has_inputs(self, cscan):
        '''
        function overridden from base class
        By definition, status = 0 if still NEED_INPUTS, -1 if NO_DATA, 1 if NEED_TO_RUN
                       qcstatus needs a value only when -1 or 0.
        You can set qcstatus to a short string that explain why it's no ready to run.
            e.g: No NIFTI

        :param cscan: object cscan define in dax.XnatUtils (see XnatUtils in dax for information)
        :return: status, qcstatus
        '''

        if XnatUtils.is_cscan_unusable(cscan):
            return -1, 'Scan unusable'

        #Check has_resource PARREC or NIFTI
        if XnatUtils.has_resource(cscan, 'NIFTI'):
            return 1, None

        LOGGER.debug('GIF Parcellation: NIFTI not found.')
        return 0, 'No NIFTI'

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
        assr_label = assessor.label()
        scan_label = assr_label.split('-x-')[3]

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   proj=proj_label,
                                   subj=subj_label,
                                   sess=sess_label,
                                   scan=scan_label,
                                   dir=jobdir,
                                   host=self.xnat_host,
                                   suffix_proc=self.suffix_proc,
                                   template=self.db_template,
                                   gif_path=self.gif)

        return [cmd]
