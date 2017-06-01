from dax import XnatUtils, spiders
import os
import shutil


# Inputs:
LABEL = '${assessor_label}'
JOB_DIR = '${temp_dir}'
B3000_FILES = '${b3000}'.split(',')
B2000_FILES = '${b2000}'.split(',')
B1500_FILES = '${b1500}'.split(',')
B500_FILES = '${b500}'.split(',')
B90_FILES = '${b90}'.split(',')
B3000_DCM = '${b3000_dcm}'.split(',')
B2000_DCM = '${b2000_dcm}'.split(',')
B1500_DCM = '${b1500_dcm}'.split(',')
B500_DCM = '${b500_dcm}'.split(',')
B90_DCM = '${b90_dcm}'.split(',')
ENV_SOURCE = '${env_source}'

# 'WIP b3000_90 SENSE', 'SWITCH DB TO YES b3000_80', 'b3000_80',
#     'b2000_vx1.3', 'b1500_vx1.3', 'b500_vx1.3', 'b90_vx1.3'


def main():
    """ Main function."""
    # Run the env_source:
    if os.path.isFile(ENV_SOURCE):
        os.system('sh {}'.format(ENV_SOURCE))

    if len(JOB_DIR) > 0 and not os.path.exists(JOB_DIR):
        os.makedirs(JOB_DIR)

    dcm_folder = XnatUtils.makedir(os.path.join(JOB_DIR, 'OsiriX'))
    fpages = list()

    # Number of acquisition:
    acquisition = 2 if len(B3000_FILES) > 1 else 1

    project, subject, session = LABEL.split('-x-')[:3]
    pdf_final = os.path.join(JOB_DIR, '%s_registration_VERDICT.pdf' % session)

    output_folder = XnatUtils.makedir(os.path.join(JOB_DIR, 'outputs'))
    osirix_folder = XnatUtils.makedir(os.path.join(output_folder, 'OsiriX'))

    for acq_ind in range(1, acquisition + 1):
        # Step 1:
        # Transform the 4D Nifti into 3D images for registration:
        print('Splitting nifti {} ...'.format(fpath))
        ori_3d = split_nifti_4D_3Ds(fpath)
        # b3000:
        ori_nii = B3000_FILES[acq_ind]
        # Convert the original NIFTI from b3000 to dicom:
        convert_nifti_2_dicom(
            ori_nii, B3000_DCM[acq_ind], B3000_DCM[acq_ind], osirix_folder,
            'b3000', label=("b3000_{}_reg".format(str(acq_ind))))

        # Other scans:
        for index, fpath in enumerate([B2000_FILES, B1500_FILES,
                                       B500_FILES, B90_FILES]):
            # Step 1:
            # Transform the 4D Nifti into 3D images for registration:
            print('Splitting nifti {} ...'.format(fpath))
            fpath_3d = split_nifti_4D_3Ds(fpath)

            if 'b3000' not in scan_info['type'].lower():
                outdir = XnatUtils.makedir(os.path.join(output_folder,
                                                        scan_info['ID']))
                # Step 2-3-4 see register nifti
                self.time_writer('Registration nifti ...')

                nii_reg = self.register_nifti(
                                    self.acquisitions[i][index],
                                    self.acquisitions[i][index-1],
                                    outdir, i)

                self.acquisitions[i][index]['reg'] = nii_reg

                # Generate DICOM version of the reg_f3d results:
                convert_nifti_2_dicom(
                    nii_reg,
                    self.acquisitions[i][index-1]['dicom'],
                    self.acquisitions[i][index]['dicom'],
                    osirix_folder,
                    scan_info['type'],
                    label=("%s_%s_reg"
                           % (scan_info['ID'],
                              scan_info['type'].replace(' ', '_'))))
            else:
                

    # Generate big niftis
    self.generate_big_nifti()

    # Make PDF
    self.make_pdf()

    # Zip the DICOMs output:
    initdir = os.getcwd()
    # Zip all the files in the directory
    zip_name = os.path.join(self.jobdir, 'outputs', 'OsiriX', 'osirix.zip')
    os.chdir(os.path.join(self.jobdir, 'outputs', 'OsiriX'))
    os.system('zip -r %s * > /dev/null' % zip_name)
    # return to the initial directory:
    os.chdir(initdir)

    # Gzip nii:
    XnatUtils.gzip_nii(os.path.join(JOB_DIR, '1', 'AMICO', '${model}'))
    if ACQ2_FILE and os.path.exists(ACQ2_FILE):
        XnatUtils.gzip_nii(os.path.join(JOB_DIR, '2', 'AMICO', '${model}'))


def make_pdf():
    """Method to make the PDF for the spider.

    :return: None
    """
    pdf_pages = dict()
    page_number = 1
    pdf_title = 'Registration Verdict - acquisition %d'
    list_slices = [3, 6, 9, 12]
    slices = {'0': list_slices,
              '1': list_slices,
              '2': list_slices,
              '3': list_slices,
              '4': list_slices}
    for i in range(1, len(self.acquisitions.keys()) + 1):
        images = list()
        labels = dict()
        sorted_list = sorted(self.acquisitions[i],
                             key=lambda k: int(k['ID']))
        for index, scan_info in enumerate(sorted_list):
            labels[str(index)] = scan_info['type']
            images.append(self.acquisitions[i][index]['reg'])

        # Saved pages:
        pdf_page = os.path.join(self.jobdir, 'registration_page_%d.pdf'
                                             % page_number)

        self.plot_images_page(pdf_page, page_number, images,
                              pdf_title % (i+1), image_labels=labels,
                              slices=slices, orient='ax')
        pdf_pages[page_number] = pdf_page
        page_number += 1
    # Merge pages:
    self.merge_pdf_pages(pdf_pages, self.pdf_final)


def register_nifti(source_info, target_info, output_folder,
                   acquisition_number):
    """Register the nifti source to the target.

    :param source_info: dictionary information of source image
    :param target_info: dictionary information of target image
    :param output_folder: path to the output folder
    :param acquisition_number: index of the acquisition
    :return: path to the 4D nifti register
    """
    # Variables:
    ala_dir = XnatUtils.makedir(os.path.join(output_folder, 'REG_ALA'))
    reg_dir = XnatUtils.makedir(os.path.join(output_folder, 'REG_RES'))
    b1tob0 = os.path.join(ala_dir, 'b1tob0.txt')
    volume_niis = {0: source_info['3D'][0]}
    # Step 2:
    # Register each scan b0 to the previous one
    # (e.g: b3000 <- b2000)
    cmd = REG_ALADIN_CMD.format(exe_path=self.reg_aladin_exe,
                                ref=target_info['3D'][0],
                                flo=source_info['3D'][0],
                                res=source_info['3D'][0],
                                aff=b1tob0,
                                args=self.args_reg_aladin)
    self.run_system_cmd(cmd)

    # Step 3:
    # Propagate the transformation to the rest of the volume:
    for index, volume in enumerate(source_info['3D']):
        if index != 0:
            vol_nii = os.path.join(reg_dir, 'volume_%s_reg.nii' % index)
            cmd = REG_ALADIN_CMD.format(exe_path=self.reg_resample_exe,
                                        ref=target_info['3D'][0],
                                        flo=volume,
                                        res=vol_nii,
                                        aff=b1tob0,
                                        args=self.args_reg_resample)
            self.run_system_cmd(cmd)
            volume_niis[index] = vol_nii

    # Step 4:
    # Put back the nifti together as one volume
    final_nii = os.path.join(output_folder, '%s_%s_%d_reg.nii'
                                            % (source_info['ID'],
                                               source_info['type'],
                                               acquisition_number))
    join_nifti_3Ds_4D(volume_niis, final_nii)
    return final_nii


def generate_big_nifti():
    """Generate big nifti with all VERDICT acquisition files."""
    for i in range(1, len(self.acquisitions.keys()) + 1):
        f_img = nib.load(self.acquisitions[i][0]['4D'])
        f_img_data = f_img.get_data()
        data = np.zeros(
                shape=(f_img_data.shape[0],
                       f_img_data.shape[1],
                       f_img_data.shape[2],
                       f_img_data.shape[3]*len(self.acquisitions[i])))

        sorted_list = sorted(self.acquisitions[i],
                             key=lambda k: int(k['ID']))
        for index, scan_info in enumerate(sorted_list):
            # Open niftis with nibabel
            f_img = nib.load(scan_info['reg'])
            f_img_data = f_img.get_data()

            for vol in range(0, f_img_data.shape[3]):
                # Draw
                vol_index = index*f_img_data.shape[3] + vol
                data[:, :, :, vol_index] = f_img_data[:, :, :, vol]

        nii_5d = nib.Nifti1Image(data, affine=f_img.affine)
        acq_dir = XnatUtils.makedir(os.path.join(self.jobdir,
                                                 'outputs',
                                                 'ACQ%d' % i))
        filename = '%s_acquisition%d.nii' % (self.xnat_session, i)
        nii_file = os.path.join(acq_dir, filename)
        nib.save(nii_5d, nii_file)
        # gzip the nifti:
        XnatUtils.gzip_file(nii_file)


def split_nifti_4D_3Ds(nifti_path):
    """Split 4D niftis into 3Ds.

    :param nifti_path: path for the nifti
    :return: list of path of new 3D niftis
    """
    # list of images:
    li_niis = list()
    # Open niftis with nibabel
    f_img = nib.load(nifti_path)
    f_img_data = f_img.get_data()
    # Draw
    if len(f_img_data.shape) < 3:
        return [nifti_path]
    elif len(f_img_data.shape) > 4:
        print('Warning: images dimension > 4. Can not be split.')
    else:
        for index in range(f_img_data.shape[3]):
            data = f_img_data[:, :, :, index]
            nii_3d = nib.Nifti1Image(data, affine=f_img.affine)
            nii_file = '%s_%d.nii' % (os.path.splitext(nifti_path)[0],
                                      index)
            nib.save(nii_3d, nii_file)
            li_niis.append(nii_file)
    return li_niis


def join_nifti_3Ds_4D(li_nii, nifti_path):
    """Join 3D niftis into 4D.

    :param li_nii: dictionary of nifti file paths in order
    :param nifti_path: nifti_path
    """
    f_img = nib.load(li_nii[0])
    f_img_data = f_img.get_data()
    data = np.zeros(shape=(f_img_data.shape[0],
                           f_img_data.shape[1],
                           f_img_data.shape[2],
                           len(li_nii.keys())))
    for index, nii in li_nii.items():
        # Open niftis with nibabel
        f_img = nib.load(nii)
        f_img_data = f_img.get_data()
        # Draw
        if len(f_img_data.shape) != 3:
            print('Warning: nifti image not 3D. Can not be join.')
        else:
            data[:, :, :, index] = f_img_data
    nii_4d = nib.Nifti1Image(data, affine=f_img.affine)
    nib.save(nii_4d, nifti_path)


def write_dicom(pixel_array, filename, ds_copy, ds_ori,
                series_number, sop_id, stype):
    """Write a dicom from a pixel_array (numpy).

    :param pixel_array: 2D numpy ndarray.
                        If pixel_array is larger than 2D, errors.
    :param filename: string name for the output file.
    :param ds_copy: pydicom object with the header that need to be copy
    :param ds_ori: original pydicom object of the pixel_array
    :param series_number: number of the series being processed
    :param sop_id: SOPInstanceUID for the DICOM
    :param stype: type for the scan
    """
    # Set the DICOM dataset
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = 'Secondary Capture Image Storage'
    file_meta.MediaStorageSOPInstanceUID = ds_ori.SOPInstanceUID
    file_meta.ImplementationClassUID = ds_ori.SOPClassUID
    ds = FileDataset(filename, {}, file_meta=file_meta, preamble="\0" * 128)

    # Copy the tag from the original DICOM
    for tag, d_obj in ds_ori.items():
        if tag != ds_ori.data_element("PixelData").tag:
            ds[tag] = d_obj

    # Other tags to set
    ds.SeriesNumber = series_number
    sop_uid = sop_id + str(datetime.datetime.now()).replace('-', '')\
                                                   .replace(':', '')\
                                                   .replace('.', '')\
                                                   .replace(' ', '')
    ds.SOPInstanceUID = sop_uid[:-1]
    ds.ProtocolName = '%s registered ' % stype
    # Set SeriesDate/ContentDate
    now = datetime.date.today()
    ds.SeriesDate = '%d%02d%02d' % (now.year, now.month, now.day)
    ds.ContentDate = '%d%02d%02d' % (now.year, now.month, now.day)
    ds.Modality = 'MR'
    ds.ConversionType = 'WSD'
    ds.StudyDescription = 'INNOVATE'
    ds.SeriesDescription = '%s_registered' % stype
    ds.SamplesperPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.SecondaryCaptureDeviceManufctur = 'Python 2.7.X'
    nb_frames = pixel_array.shape[2] * pixel_array.shape[3]
    ds.NumberOfFrames = nb_frames
    ds.PixelRepresentation = 0
    ds.HighBit = 15
    ds.BitsStored = 8
    ds.BitsAllocated = 8
    ds.SmallestImagePixelValue = pixel_array.min()
    ds.LargestImagePixelValue = pixel_array.max()
    ds.Columns = pixel_array.shape[0]
    ds.Rows = pixel_array.shape[1]

    # Fixing the sequence if the number of frames was less than the original
    # it happens if we remove the last volume for phillips data (mean)
    if ds_ori.NumberOfFrames > nb_frames:
        nb_ori_vol = ds_ori.NumberOfFrames * pixel_array.shape[3] / nb_frames
        new_seq = Sequence()
        for i in xrange(0, ds_ori.NumberOfFrames):
            if i % nb_ori_vol != 0:
                new_seq.append(ds_ori[0x5200, 0x9230][i])
        ds[0x5200, 0x9230].value = new_seq

    # Organise the array:
    size = pixel_array.shape
    array = np.zeros((size[0] * size[2] * size[3], size[1]))
    for i in range(size[2]):
        for j in range(size[3]):
            ind_down = size[0] * j + i * size[3] * size[0]
            ind_up = size[0] * (j + 1) + i * size[3] * size[0]
            array[ind_down:ind_up, :] = pixel_array[:, :, i, j]

    # Set the Image pixel array
    if array.dtype != np.uint8:
        array = array.astype(np.uint8)
    ds.PixelData = array.tostring()

    # Save the image
    ds.save_as(filename)


def convert_nifti_2_dicom(nifti_path, dcm_target, dicom_source,
                          output_folder, stype, label=None):
    """Convert 4D niftis into DICOM files.

    :param nifti_path: path to the nifti file
    :param dcm_targets: list of pydicom object from the target image for
                        header info
    :param dicom_source: one dicom file from the source for header info
    :param output_folder: folder where the DICOM files will be saved
    :param stype: type for the scan
    :param label: name for the output dicom files
    :return: None
    """
    if not os.path.isfile(nifti_path):
        raise Exception("NIFTY File %s not found." % nifti_path)
    # Load image from NIFTI
    f_img = nib.load(nifti_path)
    f_data = f_img.get_data()
    # Rotation 270
    f_data = np.rot90(f_data)
    f_data = np.rot90(f_data)
    f_data = np.rot90(f_data)
    f_data = np.flipud(f_data)

    # Scale numbers to uint8
    f_data = (f_data - f_data.min()) * 255.0 / (f_data.max() - f_data.min())

    # Load dicom headers
    if not os.path.isfile(dicom_source):
        err = "DICOM File %s not found."
        raise Exception(err % dicom_source)
    sour_obj = dicom.read_file(dicom_source)

    # Load dicom headers
    if not os.path.isfile(dcm_target):
        err = "DICOM File %s not found."
        raise Exception(err % dcm_target)
    tar_obj = dicom.read_file(dcm_target)

    # Make output_folder:
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Series Number and SOP UID
    str_ti = "%f" % time.time()
    series_number = 86532 + int(str_ti[-4:-2]) + int(str_ti[-2:])
    sop_id = sour_obj.SOPInstanceUID.split('.')
    sop_id = '.'.join(sop_id[:-1]) + '.'

    filename = os.path.join(output_folder, '%s.dcm' % (label))
    write_dicom(f_data, filename, tar_obj, sour_obj,
                series_number, sop_id, stype)


if __name__ == '__main__':
    main()
