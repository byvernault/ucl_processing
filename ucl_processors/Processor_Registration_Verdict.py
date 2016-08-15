"""Processor associated to Spider_Registration_Verdict.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Processor name: Processor_Registration_Verdict
Creation date:  2016-08-15 11:38:30.764498
Purpose:        Registration of Verdict modalities
"""

# Python packages import
import os
import logging
from dax import XnatUtils, SessionProcessor

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Registration of Verdict modalities"
__processor_name__ = "Processor_Registration_Verdict"
__modifications__ = """2016-08-15 11:38:30.764498 - Original write"""

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
HOME = os.path.expanduser("~")
DEFAULT_SPIDER_PATH = os.path.join(HOME, 'Xnat-management/ucl_processing/\
ucl_spiders/Spider_Registration_Verdict_v1_0_0.py')
DEFAULT_WALLTIME = '02:00:00'
DEFAULT_MEM = 3048
DEFAULT_VERDICT_MODALITIES = ['SWITCH DB TO YES b3000_80', 'b3000_80',
                              'b2000_vx1.3', 'b1500_vx1.3', 'b500_vx1.3',
                              'b90_vx1.3']
DEFAULT_REG_ALADIN = '/share/apps/cmic/niftypipe_deps/bin/reg_aladin'
DEFAULT_REG_RESAMPLE = '/share/apps/cmic/niftypipe_deps/bin/reg_resample'
DEFAULT_ARGS_REG_ALADIN = " -maxit 15 -ln 4 -lp 4 -interp 1"
DEFAULT_ARGS_REG_RESAMPLE = " -ln 4 -lp 4 -jl 0.1 -be 0.05 -maxit 250 \
-lncc 0 5.0 -sx 2.5"

# Format for the spider command line
SPIDER_FORMAT = """python {spider} \
-p {proj} \
-s {subj} \
-e {sess} \
-d {dir} \
--suffix "{suffix_proc}" \
--scansid {scansid} \
--regAladin {regaladin} \
--argsRegAladin "{args_reg_ala}" \
--regResample {regresample} \
--argRegResample "{args_reg_resample}" \
--openmp_core {number_cores}
"""


class Processor_Registration_Verdict(SessionProcessor):
    """Processor class for Registration_Verdict that runs on a session.

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
                 scan_modalities=DEFAULT_VERDICT_MODALITIES,
                 reg_aladin_exe=DEFAULT_REG_ALADIN,
                 args_reg_aladin=DEFAULT_ARGS_REG_ALADIN,
                 reg_resample_exe=DEFAULT_REG_RESAMPLE,
                 args_reg_resample=DEFAULT_ARGS_REG_RESAMPLE,
                 walltime=DEFAULT_WALLTIME, mem_mb=DEFAULT_MEM,
                 suffix_proc=''):
        """Entry point for Processor_Registration_Verdict Class."""
        super(Processor_Registration_Verdict,
              self).__init__(walltime, mem_mb, spider_path, version,
                             suffix_proc=suffix_proc)
        self.modalities = XnatUtils.get_input_list(scan_modalities,
                                                   DEFAULT_VERDICT_MODALITIES)
        self.reg_aladin_exe = reg_aladin_exe
        self.args_reg_aladin = args_reg_aladin
        self.reg_resample_exe = reg_resample_exe
        self.args_reg_resample = args_reg_resample

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
        # Check that there is only one scan usable with the reference type and
        # that the reference file as a NIFTI
        verdict_cscans = XnatUtils.get_good_cscans(csess, self.modalities)
        if not verdict_cscans:
            LOGGER.debug('Processor_Registration_Verdict: \
        cannot run at all, no VERDICT image found')
            return -1, 'VERDICT not found'
        for cscan in verdict_cscans:
            if not XnatUtils.has_resource(cscan, 'NIFTI'):
                LOGGER.debug('Processor_Registration_Verdict: \
        cannot run, no NIFTI found for %s scan', cscan.info()['ID'])
                return 0, "Missing NIFTI"

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

        csess = XnatUtils.CachedImageSession(assessor._intf, proj_label,
                                             subj_label, sess_label)
        source_cscans = XnatUtils.get_good_cscans(csess, self.modalities)
        sources_id = [sc.info()['ID'] for sc in source_cscans]

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   proj=proj_label,
                                   subj=subj_label,
                                   sess=sess_label,
                                   dir=jobdir,
                                   scansid=','.join(sources_id),
                                   regaladin=self.reg_aladin_exe,
                                   args_reg_ala=self.args_reg_aladin,
                                   regresample=self.reg_resample_exe,
                                   args_reg_resample=self.args_reg_resample,
                                   number_cores=self.ppn,
                                   suffix_proc=self.suffix_proc)

        return [cmd]
