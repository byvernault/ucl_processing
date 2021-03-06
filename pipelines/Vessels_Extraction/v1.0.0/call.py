import os
import glob
from dax import spiders, XnatUtils


JOB_DIR = '${temp_dir}'
MIN = '${min}'
MAX = '${max}'
CT = '${ct}'
VESSELS_FILES = '${vessels}'
NIFTYPIPE_EXE = '${exe}'
OPENMP_CORE = '${omp}'
WORKING_DIR = '${working_dir}'
EXE_CMD = '''{exe_path} \
--vessels {vessels} \
--output {output} \
--no_qsub \
--remove_tmp \
--n_procs 1 \
{ct} \
{min} \
{max} \
{omp} \
{wdir}'''
OMP = '''--openmp_core {number_core}'''


def main():
    """ Main function."""
    _working_dir = None
    if WORKING_DIR != 'None':
        _working_dir = os.path.join(WORKING_DIR, '${assessor_label}')
        if not os.path.exists(_working_dir):
            os.makedirs(_working_dir)

    if XnatUtils.executable_exists(NIFTYPIPE_EXE):
        _omp = '' if OPENMP_CORE == 'None' \
               else OMP.format(number_core=OPENMP_CORE)
        _wd = '' if _working_dir is not None \
              else '--working_dir {0}'.format(_working_dir)
        _min = '' if MIN == 'None' else '--min {0}'.format(MIN)
        _max = '' if MAX == 'None' else '--max {0}'.format(MAX)

        cmd = EXE_CMD.format(exe_path=NIFTYPIPE_EXE,
                             vessels=" ".join(VESSELS_FILES.split(',')),
                             output=JOB_DIR,
                             omp=_omp,
                             wdir=_wd,
                             min=_min,
                             max=_max,
                             ct='--ct' if CT == '1' else '')
        os.system(cmd)
        make_pdf()
    else:
        raise Exception("Error: executable %s not found" % (NIFTYPIPE_EXE))


def make_pdf():
    """Method to make the PDF for the spider.

    :return: None
    """
    # PDF pages:
    pdf = os.path.join(JOB_DIR, 'vessels_extraction.pdf')

    # Images outputs:
    images = VESSELS_FILES.split(',')
    images.extend(glob.glob(os.path.join(JOB_DIR, '*_extracted.nii.gz')))

    labels = dict()
    for ind in range(len(images)):
        if ind == 0 or ind == 1:
            labels[str(ind)] = 'Vessel {0}'.format(str(ind))
        else:
            labels[str(ind)] = 'Vessel Extracted {0}'.format(str(ind))

    spiders.plot_images(pdf, 1, images, 'Vessel Extraction Pipeline',
                        image_labels=labels)


if __name__ == '__main__':
    main()
