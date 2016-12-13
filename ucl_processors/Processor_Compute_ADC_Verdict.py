"""Processor associated to Spider_Compute_ADC_Verdict.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Processor name: Processor_Compute_ADC_Verdict
Creation date:  2016-08-18 13:14:30.689905
Purpose:        Generate ADC map from Verdict registered scans
"""

# Python packages import
import os
import logging
from dax import XnatUtils, SessionProcessor

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
ucl_spiders/', 'Spider_Compute_ADC_Verdict_v1_0_0.py')
DEFAULT_WALLTIME = '00:30:00'
DEFAULT_MEM = 6048
DEFAULT_PROCTYPE = 'Registration_Verdict_v1'
DEFAULT_MATLAB_CODE = os.path.join(HOME, 'Code', 'matlab', 'VERDICT',
                                   'ADC_MAP')
DEFAULT_CAMINO = os.path.join(HOME, 'Code', 'caminoLaura')
DEFAULT_VERDICT_MODALITIES = [
    'WIP b3000_90 SENSE', 'SWITCH DB TO YES b3000_80', 'b3000_80',
    'b2000_vx1.3', 'b1500_vx1.3', 'b500_vx1.3', 'b90_vx1.3']
DEFAULT_SCHEME_FILE = os.path.join(HOME, 'Code', 'matlab', 'VERDICT',
                                   'Noptimised', 'NOptimisedADC_IN.scheme')

# Format for the spider command line
SPIDER_FORMAT = """python {spider} \
-p {proj} \
-s {subj} \
-e {sess} \
-d {dir} \
--suffix "{suffix_proc}" \
--proctype {proctype} \
--nbAcq {nb_acq} \
--mc {mc} \
--camino {camino} \
--scheme {scheme} \
"""


class Processor_Compute_ADC_Verdict(SessionProcessor):
    """Processor class for Verdict that runs on a session.

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider

    :param suffix: suffix to the spider
    """

    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 scan_modalities=DEFAULT_VERDICT_MODALITIES,
                 proctype=DEFAULT_PROCTYPE, matlab_code=DEFAULT_MATLAB_CODE,
                 camino=DEFAULT_CAMINO, scheme_file=DEFAULT_SCHEME_FILE,
                 walltime=DEFAULT_WALLTIME,
                 mem_mb=DEFAULT_MEM,  suffix_proc=''):
        """Entry point for Processor_Verdict Class."""
        super(Processor_Compute_ADC_Verdict,
              self).__init__(walltime, mem_mb, spider_path, version,
                             suffix_proc=suffix_proc)
        self.modalities = XnatUtils.get_input_list(scan_modalities,
                                                   DEFAULT_VERDICT_MODALITIES)
        self.proctype = proctype
        self.mc = matlab_code
        self.camino = camino
        self.scheme_file = scheme_file

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
        verdict_cscans = XnatUtils.get_good_cscans(csess, self.modalities)
        if not verdict_cscans:
            LOGGER.debug('Processor_Registration_Verdict: \
        cannot run at all, no VERDICT image found')
            return -1, 'VERDICT not found'

        verdict_cassrs = list()
        for cassr in csess.assessors():
            if XnatUtils.is_cassessor_good_type(cassr, [self.proctype]):
                verdict_cassrs.append(cassr)
        if not verdict_cassrs:
            LOGGER.debug('Processor_Compute_ADC_Verdict: \
        cannot run, no good QA Registration VERDICT found')
            return 0, 'Registration missing'

        cassr = verdict_cassrs[0]
        LOGGER.debug('Processor_Compute_ADC_Verdict: \
good registration assessor found: %s', cassr.info()['label'])

        if not XnatUtils.has_resource(cassr, 'ACQ1'):
            LOGGER.debug('Processor_Compute_ADC_Verdict: \
cannot run, no ACQ resource found for %s assessor',
                         cassr.info()['label'])
            return 0, "Missing ACQ#"

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

        nb_acq = 1
        csess = XnatUtils.CachedImageSession(assessor._intf, proj_label,
                                             subj_label, sess_label)
        reg_verdict = ''
        for cassr in csess.assessors():
            if XnatUtils.is_cassessor_good_type(cassr, [self.proctype]):
                reg_verdict = cassr

        if XnatUtils.has_resource(reg_verdict, 'ACQ2'):
            nb_acq = 2

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   proj=proj_label,
                                   subj=subj_label,
                                   sess=sess_label,
                                   dir=jobdir,
                                   nb_acq=nb_acq,
                                   proctype=self.proctype,
                                   mc=self.mc,
                                   camino=self.camino,
                                   scheme=self.scheme_file,
                                   suffix_proc=self.suffix_proc)

        return [cmd]
