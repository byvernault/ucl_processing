"""Processor associated to Spider_BaMoS.

Author:         Benjamin Yvernault
contact:        b.yvernault@ucl.ac.uk
Processor name: Processor_BaMoS
Creation date:  2016-12-21 15:42:01
Purpose:        Generate Assessor for BaMoS Spider
"""

# Python packages import
import os
import logging
from dax import XnatUtils, SessionProcessor

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Generate Assessor for BaMoS Spider."
__processor_name__ = "Processor_BaMoS"
__modifications__ = """2016-12-21 15:42:01 - Original write"""

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
HOME = os.path.expanduser("~")
DEFAULT_SPIDER_PATH = os.path.join(HOME, 'Xnat-management', 'ucl_processing',
                                   'ucl_spiders', 'Spider_BaMoS_v1_0_1.py')
DEFAULT_WALLTIME = '10:00:00'
DEFAULT_MEM = 8048
DEFAULT_BAMOS_SCRIPT = os.path.join(HOME, 'Code', 'scripts', 'bash',
                                    'BaMoSGenericDax.sh')
DEFAULT_GIF_PROCTYPE = 'GIF_Parcellation_v2'
DEFAULT_REG_FOLDER = '/share/apps/cmic/niftypipe_deps/bin/'
DEFAULT_SEG_FOLDER = '/share/apps/cmic/niftypipe_deps/bin/'
DEFAULT_T1 = ['T1', 'MPRAGE']
DEFAULT_FLAIR = ['FLAIR']
DEFAULT_T2 = ['T2']
DEFAULT_SEGBIASM = '/home/csudre/NiftySeg_0.9.4/build_comic2/seg-apps/\
Seg_BiASM'
DEFAULT_SEGANA = '/home/csudre/NiftySeg_0.9.4/build_comic2/seg-apps/\
Seg_Analysis'
DEFAULT_GMATRIX = '/cluster/project0/SegBiASM/DataToTryBaMoS/GMatrix4_Low3.txt'
DEFAULT_RULE = '/cluster/project0/SegBiASM/DataToTryBaMoS/GenericRule_CSF.txt'

# Format for the spider command line
SPIDER_FORMAT = """python {spider} \
-p {proj} \
-s {subj} \
-e {sess} \
-d {dir} \
--suffix "{suffix_proc}" \
--t1 {t1} \
--flair {flair} \
--gif {gif} \
--bamos {bamos} \
--regfolder {reg} \
--segfolder {seg} \
--seg_biasm {biasm} \
--seg_analysis {ana} \
--gmatrix {gmatrix} \
--rule {rule}
"""


class Processor_BaMoS(SessionProcessor):
    """Processor class for Verdict that runs on a session.

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider

    :param suffix: suffix to the spider
    """

    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None,
                 t1=DEFAULT_T1, flair=DEFAULT_FLAIR, t2=DEFAULT_T2,
                 gif_proctype=DEFAULT_GIF_PROCTYPE,
                 bamos_script=DEFAULT_BAMOS_SCRIPT,
                 reg_folder=DEFAULT_REG_FOLDER, seg_folder=DEFAULT_SEG_FOLDER,
                 seg_biasm=DEFAULT_SEGBIASM, seg_analysis=DEFAULT_SEGANA,
                 gmatrix=DEFAULT_GMATRIX, rule=DEFAULT_RULE,
                 walltime=DEFAULT_WALLTIME, mem_mb=DEFAULT_MEM,
                 suffix_proc=''):
        """Entry point for Processor_Verdict Class."""
        super(Processor_BaMoS,
              self).__init__(walltime, mem_mb, spider_path, version,
                             suffix_proc=suffix_proc)
        self.gif = gif_proctype
        self.t1 = XnatUtils.get_input_list(t1, DEFAULT_T1)
        self.flair = XnatUtils.get_input_list(flair, DEFAULT_FLAIR)
        self.t2 = XnatUtils.get_input_list(t2, DEFAULT_T2)
        self.bamos = bamos_script
        self.regfolder = reg_folder
        self.segfolder = seg_folder
        self.seg_biasm = seg_biasm
        self.seg_analysis = seg_analysis
        self.gmatrix = gmatrix
        self.rule = rule

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
        t1_scans = XnatUtils.get_good_cscans(csess, self.t1)
        if not t1_scans:
            LOGGER.debug('Processor_BaMoS: \
        cannot run at all, no T1 image found')
            return -1, 'T1 not found'

        flair_scans = XnatUtils.get_good_cscans(csess, self.flair)
        if not flair_scans:
            LOGGER.debug('Processor_BaMoS: \
        cannot run at all, no FLAIR image found')
            return -1, 'FLAIR not found'

        gif_proctype = list()
        for cassr in csess.assessors():
            if XnatUtils.is_cassessor_good_type(cassr, [self.gif]):
                gif_proctype.append(cassr)
        if not gif_proctype:
            LOGGER.debug('Processor_BaMoS: \
        cannot run, no good QA %s found' % self.gif)
            return 0, '%s missing' % self.gif

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

        t1_scans = XnatUtils.get_good_cscans(csess, self.t1)
        flair_scans = XnatUtils.get_good_cscans(csess, self.flair)
        t2_scans = XnatUtils.get_good_cscans(csess, self.t2)

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   proj=proj_label,
                                   subj=subj_label,
                                   sess=sess_label,
                                   dir=jobdir,
                                   t1=t1_scans[0].info()['ID'],
                                   flair=flair_scans[0].info()['ID'],
                                   gif=self.gif,
                                   bamos=self.bamos,
                                   reg=self.regfolder,
                                   seg=self.segfolder,
                                   biasm=self.seg_biasm,
                                   ana=self.seg_analysis,
                                   rule=self.rule,
                                   gmatrix=self.gmatrix,
                                   suffix_proc=self.suffix_proc)

        # Add the T2 if found
        if t2_scans:
            cmd = '%s --t2 %s' % t2_scans[0].info()['ID']

        return [cmd]
