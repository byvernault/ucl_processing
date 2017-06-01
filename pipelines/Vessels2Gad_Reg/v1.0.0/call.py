import os
import glob
from dax import spiders, XnatUtils


JOB_DIR = '${temp_dir}'
GAD_FILE = '${gad}'
VESSEL_FILE = '${vessel}'
NIFTYPIPE_EXE = '${exe}'
OPENMP_CORE = '${omp}'
WORKING_DIR = '${working_dir}'
EXE_CMD = '''{exe_path} \
--ref {gad} \
--vessels {vessels} \
--output {output} \
--no_qsub \
--remove_tmp \
--n_procs 1 \
{omp} \
{wdir}'''
OMP = '''--openmp_core {number_core}'''
WDIR = '''--working_dir '{working_dir}' '''


def main():
    """ Main function."""
    if WORKING_DIR != 'None' and not os.path.exists(WORKING_DIR):
        os.makedirs(WORKING_DIR)

    if XnatUtils.executable_exists(NIFTYPIPE_EXE):
        _omp = ''
        _wd = ''
        if OPENMP_CORE is not None and OPENMP_CORE != 'None':
            _omp = OMP.format(number_core=OPENMP_CORE)
        if WORKING_DIR is not None and WORKING_DIR != 'None':
            _wd = WDIR.format(working_dir=WORKING_DIR)

        cmd = EXE_CMD.format(exe_path=NIFTYPIPE_EXE,
                             gad=GAD_FILE,
                             vessels=VESSEL_FILE,
                             output=JOB_DIR,
                             omp=_omp,
                             wdir=_wd)
        os.system(cmd)
        make_pdf()
    else:
        raise Exception("Error: executable %s not found" % (NIFTYPIPE_EXE))


def make_pdf():
    """Method to make the PDF for the spider.

    :return: None
    """
    # PDF pages:
    pdf = os.path.join(JOB_DIR, 'vessels2gad.pdf')

    # Images outputs:
    images = [GAD_FILE]
    images.extend(glob.glob(os.path.join(JOB_DIR, '*.nii.gz')))

    labels = dict()
    for ind in range(len(images)):
        if ind == 0:
            labels[str(ind)] = 'Reference'
        else:
            labels[str(ind)] = 'Vessel {0}'.format(str(ind))

    spiders.plot_images(pdf, 1, images, 'Vessel to Gad Registration Pipeline',
                        image_labels=labels)


if __name__ == '__main__':
    main()
