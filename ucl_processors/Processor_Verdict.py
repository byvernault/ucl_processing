"""Processor associated to Spider_Verdict.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Processor name: Processor_Verdict
Creation date:  2016-08-18 13:14:30.689905
Purpose:        Generate Verdict Map from all Verdict scans
"""

# Python packages import
import os
import logging
from dax import XnatUtils, SessionProcessor

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Generate Verdict Map from all Verdict scans"
__processor_name__ = "Processor_Verdict"
__modifications__ = """2016-08-18 13:14:30.689905 - Original write"""

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
HOME = os.path.expanduser("~")
DEFAULT_SPIDER_PATH = os.path.join(HOME, 'Xnat-management/ucl_processing/\
ucl_spiders/', 'Spider_Verdict_v1_0_0.py')
DEFAULT_WALLTIME = '02:00:00'
DEFAULT_MEM = 3048
DEFAULT_PROCTYPE = ['Registration_Verdict_v1']
DEFAULT_MATLAB_CODE = os.path.join(HOME, 'Code', 'matlab')
DEFAULT_AMICO = os.path.join(HOME, 'Code', 'AMICO', 'matlab')
DEFAULT_CAMINO = os.path.join(HOME, 'Code', 'caminoLaura', 'bin')
DEFAULT_SPAMS = os.path.join(HOME, 'Code', 'spams-matlab')

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
--amico {amico} \
--camino {camino} \
--spams {spams} \
"""


class Processor_Verdict(SessionProcessor):
    """Processor class for Verdict that runs on a session.

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider

    :param suffix: suffix to the spider
    """

    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 proctype=DEFAULT_PROCTYPE, matlab_code=DEFAULT_MATLAB_CODE,
                 amico=DEFAULT_AMICO, camino=DEFAULT_CAMINO,
                 spams=DEFAULT_SPAMS, walltime=DEFAULT_WALLTIME,
                 mem_mb=DEFAULT_MEM,  suffix_proc=''):
        """Entry point for Processor_Verdict Class."""
        super(Processor_Verdict,
              self).__init__(walltime, mem_mb, spider_path, version,
                             suffix_proc=suffix_proc)
        self.nb_acq = 1
        self.proctype = XnatUtils.get_input_list(proctype, DEFAULT_PROCTYPE)
        self.mc = matlab_code
        self.amico = amico
        self.camino = camino
        self.spams = spams

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
        verdict_csassr = XnatUtils.get_good_cassr(csess, self.proctype)
        if not verdict_csassr:
            LOGGER.debug('Processor_Verdict: \
        cannot run at all, no Registration VERDICT found')
            return -1, 'VERDICT not found'
        for cassr in verdict_csassr:
            if not XnatUtils.has_resource(cassr, 'ACQ1'):
                LOGGER.debug('Processor_Registration_Verdict: \
cannot run, no ACQ resource found for %s assessor',
                             cassr.info()['label'])
                return 0, "Missing ACQ#"
            if XnatUtils.has_resource(cassr, 'ACQ2'):
                self.nb_acq = 2

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

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   proj=proj_label,
                                   subj=subj_label,
                                   sess=sess_label,
                                   dir=jobdir,
                                   nb_acq=self.nb_acq,
                                   proctype=self.proctype,
                                   mc=self.mc,
                                   amico=self.amico,
                                   camino=self.camino,
                                   spams=self.spams)

        return [cmd]
