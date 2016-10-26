"""Processor associated to Spider_Diffusion_Model_Fitting.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Processor name: Processor_Diffusion_Model_Fitting
Creation date:  2016-05-11 11:58:34.817004
Purpose:        Processor for Diffusion Model Fitting with pre-processing
                steps.
"""

# Python packages import
import os
import logging
from dax import XnatUtils, SessionProcessor

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Processor for Diffusion Model Fitting with \
pre-processing steps."
__processor_name__ = "Processor_Diffusion_Model_Fitting"
__modifications__ = """2016-05-11 11:58:34.817004 - Original write"""

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
HOME = os.path.expanduser("~")
DEFAULT_SPIDER_PATH = os.path.join(HOME, 'Xnat-management/ucl_processing/\
ucl_spiders/Spider_Diffusion_Model_Fitting_v1_0_0.py')
DEFAULT_WALLTIME = '10:00:00'
DEFAULT_MEM = 4096
DEFAULT_PPN = 4
DEFAULT_DTI_TYPES = ['DTI', 'DWI']
DEFAULT_T1_TYPES = ['T1', 'MPRAGE']
DEFAULT_GIF_TYPE = ['GIF_Parcellation_v1']
DEFAULT_DTI_PATH = os.path.join(HOME,
                                'anaconda/bin/perform_dti_processing.py')
DEFAULT_DTI_ARGS = "--rot 34.56 --etd 2.46 --ped -y"
DEFAULT_WORKING_DIR = '/scratch0/dax/'

# Format for the spider command line
SPIDER_FORMAT = """python {spider} \
-p {proj} \
-s {subj} \
-e {sess} \
-d {dir} \
--suffix "{suffix_proc}" \
--dti {dti} \
--gif {gif} \
--exe {exe} \
--dtiargs {exe_args} \
--openmp_core {number_core} \
--working_dir "{working_dir}"
"""


class Processor_Diffusion_Model_Fitting(SessionProcessor):
    """Processor class for Diffusion_Model_Fitting that runs on a session.

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider
    :param dtitypes: dti types on XNAT
    :param t1types: t1 types on XNAT given to GIF
    :param giftype: GIF proctype that ran on XNAT
    :param exe: path to perform_dti_processing executable
    :param exe_args: arguments for perform_dti_processing (rot/etd/ped)
    :param ppn: number of core used by the jobs
    :param suffix: suffix to the spider
    :param working_dir: working directory for spider temp files
    """

    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 walltime=DEFAULT_WALLTIME, mem_mb=DEFAULT_MEM,
                 ppn=DEFAULT_PPN, dtitypes=DEFAULT_DTI_TYPES,
                 t1types=DEFAULT_T1_TYPES, giftypes=DEFAULT_GIF_TYPE,
                 exe=DEFAULT_DTI_PATH, exe_args=DEFAULT_DTI_ARGS,
                 suffix_proc='', working_dir=DEFAULT_WORKING_DIR):
        """Entry point for Processor_Diffusion_Model_Fitting Class."""
        super(Processor_Diffusion_Model_Fitting,
              self).__init__(walltime, mem_mb, spider_path, version, ppn=ppn,
                             suffix_proc=suffix_proc)
        # DTI types:
        self.dtitypes = XnatUtils.get_input_list(dtitypes, DEFAULT_DTI_TYPES)
        # T1 types:
        self.t1types = XnatUtils.get_input_list(t1types, DEFAULT_T1_TYPES)
        # GIF proctype:
        self.giftypes = XnatUtils.get_input_list(giftypes, DEFAULT_GIF_TYPE)
        # Exe gif + args:
        self.exe = exe
        self.exe_args = exe_args
        # working dir:
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
        # DTI:
        dti_cscans = XnatUtils.get_good_cscans(csess, self.dtitypes)
        if not dti_cscans:
            LOGGER.debug('Processor_Diffusion_Model_Fitting: \
cannot run at all, no DTI images found')
            return -1, 'DTI not found'
        for dti_cscan in dti_cscans:
            if not XnatUtils.has_resource(dti_cscan, 'NIFTI'):
                LOGGER.debug('Processor_Diffusion_Model_Fitting: cannot run, \
no NIFTI found for DTI image %s' % dti_cscan.info()['ID'])
                return 0, "no NIFTI for %s" % dti_cscan.info()['ID']
        # T1:
        t1_cscans = XnatUtils.get_good_cscans(csess, self.t1types)
        if not t1_cscans:
            LOGGER.debug('Processor_Diffusion_Model_Fitting: \
cannot run at all, no T1 images found')
            return -1, 'T1 not found'
        # GIF:
        gif_cassrs = XnatUtils.get_good_cassr(csess, self.giftypes)
        if not gif_cassrs:
            LOGGER.debug('Processor_Diffusion_Model_Fitting: \
cannot run, no good GIF Parcellation found')
            return 0, 'good GIF not found'
        if len(gif_cassrs) > 1:
            LOGGER.debug('Processor_Diffusion_Model_Fitting: \
cannot run, no good GIF Parcellation found')
            return 0, 'too many good GIF'

        return 1, None

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

        csess = XnatUtils.CachedImageSession(assessor._intf, proj_label,
                                             subj_label, sess_label)
        dti_cscans = XnatUtils.get_good_cscans(csess, self.dtitypes)
        dtis = ','.join([cscan.info()['ID'] for cscan in dti_cscans])
        gif_cassrs = XnatUtils.get_good_cassr(csess, self.giftypes)
        working_dir = os.path.join(self.working_dir, assr_label)

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   proj=proj_label,
                                   subj=subj_label,
                                   sess=sess_label,
                                   dir=jobdir,
                                   suffix_proc=self.suffix_proc,
                                   dti=dtis,
                                   gif=gif_cassrs[0].info()['label'],
                                   exe=self.exe,
                                   gif_path=self.exe_args,
                                   number_core=self.ppn,
                                   working_dir=working_dir)

        return [cmd]
