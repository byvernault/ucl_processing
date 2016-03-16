'''
    Author:         Benjamin Yvernault
    contact:        b.yvernault@ucl.ac.uk
    Processor name: Processor_Reg_ADC_2_T2
    Creation date:  2016-03-16 10:17:34.887155
    Purpose:        Processor for registering ADC scan to T2
'''

__author__ = "Benjamin Yvernault"
__email__ = "b.yvernault@ucl.ac.uk"
__purpose__ = "Processor for registering ADC scan to T2"
__processor_name__ = "Processor_Reg_ADC_2_T2"
__modifications__ = "2016-03-16 10:17:34.887155 - Original write"

# Python packages import
import logging
from dax import XnatUtils, SessionProcessor

# set-up logger for printing statements
LOGGER = logging.getLogger('dax')

# Default values for arguments:
# EDIT PARAMETERS FOR YOUR SPIDER CASE (SPIDER_PATH, WALLTIME, etc...)
DEFAULT_SPIDER_PATH = '/home/byvernau/Xnat-management/ucl_processing/ucl_spiders/Spider_Reg_ADC_2_T2_v1_0_0.py'
DEFAULT_WALLTIME = '02:00:00'
DEFAULT_MEM = 2048
DEFAULT_PPN = 1
DEFAULT_T2 = ['T2 Axial']
DEFAULT_ADC = ['ADC']
DEFAULT_REG_ALADIN = '/share/apps/cmic/niftypipe_deps/bin/reg_aladin'
DEFAULT_REG_F3D = '/share/apps/cmic/niftypipe_deps/bin/reg_f3d'
DEFAULT_ARGS_REG_ALADIN = " -maxit 15 -ln 4 -lp 4 -interp 1"
DEFAULT_ARGS_REG_F3D = " -ln 4 -lp 4 -jl 0.1 -be 0.05 -maxit 250 -lncc 0 5.0 -sx 2.5"

# Format for the spider command line
SPIDER_FORMAT = '''python {spider} -p {proj} -s {subj} -e {sess} -d {dir} --suffix "{suffix_proc}" --adc {adc} --t2 {t2} --regAladin {regaladin} --argsRegAladin "{args_reg_ala}" --regf3d {regf3d} --argRegf3d "{args_reg_f3d}" --openmp_core {number_cores}'''

class Processor_Reg_ADC_2_T2(SessionProcessor):
    '''
    Processor class for Reg_ADC_2_T2 that runs on a session

    :param spider_path: spider path on the system
    :param version: version of the spider
    :param walltime: walltime required by the spider
    :param mem_mb: memory in Mb required by the spider
    :param ppn: number of cores

    :param suffix: suffix to the spider
    '''
    def __init__(self, spider_path=DEFAULT_SPIDER_PATH, version=None, adc=DEFAULT_ADC, t2=DEFAULT_T2,
                 reg_aladin_exe=DEFAULT_REG_ALADIN, args_reg_aladin=DEFAULT_ARGS_REG_ALADIN,
                 reg_f3d_exe=DEFAULT_REG_F3D, args_reg_f3d=DEFAULT_ARGS_REG_F3D,
                 walltime=DEFAULT_WALLTIME, mem_mb=DEFAULT_MEM, ppn=DEFAULT_PPN,
                 suffix_proc=''):
        super(Processor_Reg_ADC_2_T2, self).__init__(walltime, mem_mb, spider_path, version,
                                                     ppn=ppn, suffix_proc=suffix_proc)
        self.t2 = t2
        self.adc = adc
        self.reg_aladin_exe = reg_aladin_exe
        self.args_reg_aladin = args_reg_aladin
        self.reg_f3d_exe = reg_f3d_exe
        self.args_reg_f3d = args_reg_f3d

    def has_inputs(self, csess):
        '''
        function overridden from base class
        By definition, status = 0 if still NEED_INPUTS, -1 if NO_DATA, 1 if NEED_TO_RUN
                       qcstatus needs a value only when -1 or 0.
        You can set qcstatus to a short string that explain why it's no ready to run.
            e.g: No NIFTI

        :param csess: object csess define in dax.XnatUtils (see XnatUtils in dax for information)
        :return: status, qcstatus
        '''
        # Check that there is only one scan usable with the reference type and that the reference file as a NIFTI
        target_cscans = XnatUtils.get_good_cscans(csess, self.t2)
        if not target_cscans:
            LOGGER.debug('Processor_Registration2Ref: cannot run at all, no T2 image found')
            return -1, 'T2 not found'
        if len(target_cscans) > 1:
            LOGGER.debug('Processor_Registration2Ref: cannot run at all, too many T2 images found')
            return 0, 'Too many T2 scans'
        if not XnatUtils.has_resource(target_cscans[0], 'NIFTI'):
            LOGGER.debug('Processor_Registration2Ref: cannot run, no NIFTI for T2 image')
            return 0, "no T2's NIFTI"

        source_cscans = XnatUtils.get_good_cscans(csess, self.adc)
        if not source_cscans:
            LOGGER.debug('Processor_Registration2Ref: cannot run at all, no ADC image found')
            return -1, 'ADC not found'
        if len(source_cscans) > 1:
            LOGGER.debug('Processor_Registration2Ref: cannot run at all, too many ADC images found')
            return 0, 'Too many ADC scans'
        if not XnatUtils.has_resource(source_cscans[0], 'NIFTI'):
            LOGGER.debug('Processor_Registration2Ref: cannot run, no NIFTI for ADC image')
            return 0, "no ADC's NIFTI"

        return 1, None

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

        csess = XnatUtils.CachedImageSession(assessor._intf, proj_label, subj_label, sess_label)
        target_cscans = XnatUtils.get_good_cscans(csess, self.t2)
        target_id = target_cscans[0].info()['ID']
        source_cscan = XnatUtils.get_good_cscans(csess, self.adc)
        source_id = source_cscan[0].info()['ID']

        cmd = SPIDER_FORMAT.format(spider=self.spider_path,
                                   proj=proj_label,
                                   subj=subj_label,
                                   sess=sess_label,
                                   dir=jobdir,
                                   t2=target_id,
                                   adc=source_id,
                                   regaladin=self.reg_aladin_exe,
                                   args_reg_ala=self.args_reg_aladin,
                                   regf3d=self.reg_f3d_exe,
                                   args_reg_f3d=self.args_reg_f3d,
                                   number_cores=self.ppn,
                                   suffix_proc=self.suffix_proc)

        return [cmd]
