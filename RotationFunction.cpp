//
// Created by cxfeng on 19-2-21.
//

#include "RotationFunction.h"
#include "RotateData.h"

RotationFunction::RotationFunction(){

}

RotationFunction::~RotationFunction(){

}

#if 1
void RotationFunction::inititializationData(){

    pelDir = "/home/cxfeng/Desktop/My_Work/BinaryData/Binary2/pelVol.bin";
    leftDir = "/home/cxfeng/Desktop/My_Work/BinaryData/Binary2/leftVol.bin";
    rightDir = "/home/cxfeng/Desktop/My_Work/BinaryData/Binary2/rightVol.bin";
    dicomDir = "/home/cxfeng/Desktop/My_Work/DicomData/12";

    dicomDimension = new int [3]{512, 512, 600};
    spacings = new double [3]{0.96875, 0.96875, 1};
    pngSaveName = "pngSaveName";
    independentComponents = true;

    pelRotateAngle = new double [3]{-8.19445, 0.0733765, 3.75134};
    pelRotateCenter = new double [3]{155, 211.188, 83};
    leftRotateAngle = new double [3]{10.5699, 0.484976, 18.9191};
    leftRotateCenter = new double [3]{340.088, 251.34, 157.613};
    rightRotateAngle = new double [3]{10.8356, -0.788843, -7.88475};
    rightRotateCenter = new double [3]{163.947, 255.421, 166.661};

    cameraPosition = new double [3]{};
    cameraViewUp = new double [3]{};
    cameraFocalPoint = new double [3]{};
    cameraClippingRange = new double [2]{};


    vtkUcharArray_Pel =
            vtkSmartPointer<vtkUnsignedCharArray>::New();
    vtkUcharArray_Left =
            vtkSmartPointer<vtkUnsignedCharArray>::New();
    vtkUcharArray_Right =
            vtkSmartPointer<vtkUnsignedCharArray>::New();

    // 4 imagedata
    pelImageData =
            vtkSmartPointer<vtkImageData>::New();
    leftImageData =
            vtkSmartPointer<vtkImageData>::New();
    rightImageData =
            vtkSmartPointer<vtkImageData>::New();
    dicomImageData =
            vtkSmartPointer<vtkImageData>::New();

    // 4 mapper
    pelVolumeMapper =
            vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
    leftVolumeMapper =
            vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
    rightVolumeMapper =
            vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
    dicomVolumeMapper =
            vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
    // 4 volume
    pelVolume =
            vtkSmartPointer<vtkVolume>::New();
    leftVolume =
            vtkSmartPointer<vtkVolume>::New();
    rightVolume =
            vtkSmartPointer<vtkVolume>::New();
    dicomVolume =
            vtkSmartPointer<vtkVolume>::New();

    renderWindow =
            vtkSmartPointer<vtkRenderWindow>::New();
    interactorStyle =
            vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    renderWindowInteractor =
            vtkSmartPointer<vtkRenderWindowInteractor>::New();
    volumeProperty =
            vtkSmartPointer<vtkVolumeProperty>::New();
    volumeOpacity =
            vtkSmartPointer<vtkPiecewiseFunction>::New();
    color =
            vtkSmartPointer<vtkColorTransferFunction>::New();
    Renderer =
            vtkSmartPointer<vtkRenderer>::New();

    oriAxesActor =
            vtkSmartPointer<vtkAxesActor>::New();
    widget =
            vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    camera =
            vtkSmartPointer<vtkCamera>::New();
}
#endif
void RotationFunction::readDicomSeries(std::string dataDir, vtkImageData * ptr){
    vtkSmartPointer<vtkDICOMImageReader> dicomReader =
            vtkSmartPointer<vtkDICOMImageReader>::New();
    dicomReader->SetDirectoryName(dataDir.c_str());
    dicomReader->Update();

    ptr = dicomReader->GetOutput();
}

void RotationFunction::readMaskBinary(std::string dataDir, vtkUnsignedCharArray* ptr1, vtkImageData * ptr2){

    unsigned char *U_charArray_P = new unsigned char[dicomSize];
    char *binArray_P = new char[dicomSize];

    ifstream binFile(dataDir.c_str(), ios::in | ios::binary);

    binFile.read(binArray_P, dicomSize);
    binFile.close();
    U_charArray_P = reinterpret_cast<unsigned char *>(binArray_P);
    ptr1->SetArray(U_charArray_P, dicomSize, 1);

    //vtkImageData * tempImageData = vtkImageData::New();
    ptr2->GetPointData()->SetScalars(ptr1);
    ptr2->GetScalarSize();

    ptr2->SetDimensions(dicomDimension[0], dicomDimension[1],dicomDimension[2]);
    ptr2->SetSpacing(spacings[0], spacings[1], spacings[2]);
    ptr2->Modified();

}


void RotationFunction::maskAndDicomAreMultipliedToDivideWhenYAxisMirror(vtkImageData* imageData1, vtkImageData* imageData2){
    //If the Y axis is mirror symmetrical
    unsigned char *ptr = (unsigned char *) imageData1->GetScalarPointer();
    for (int z = 0; z < dicomDimension[2]; ++z){
        for (int y = 0; y < dicomDimension[1]; ++y){
            for (int x = 0; x < dicomDimension[0]; ++x){

                double *pixel = static_cast<double *>(imageData2->GetScalarPointer(x, y, z));
                *pixel = (*pixel) * (*(ptr + x + ((dicomDimension[1]-y-1)*dicomDimension[0]) + (z*dicomDimension[1]*dicomDimension[0])));

            }
        }
    }
}

void RotationFunction::maskAndDicomAreMultipliedToDivide(vtkImageData* imageData1, vtkImageData* imageData2){

    for (int z = 0; z < dicomDimension[2]; ++z){
        for (int y = 0; y < dicomDimension[1]; ++y){
            for (int x = 0; x < dicomDimension[0]; ++x){
                double *pixel = static_cast<double *>(imageData2->GetScalarPointer(x, y, z));
                unsigned char *ptr = (unsigned char *) imageData1->GetScalarPointer(x, y, z);
                *pixel = (*pixel) * (*ptr);
            }
        }
    }
}

void RotationFunction::rotateVolume(vtkVolume* volume, double rotateCenter[3], double rotateAngle[3]){
    vtkSmartPointer<vtkTransform> transform =
            vtkSmartPointer<vtkTransform>::New();
    transform->PostMultiply();
    transform->Translate(-rotateCenter[0], -rotateCenter[1], -rotateCenter[2]);
    transform->RotateX(rotateAngle[0]);
    transform->RotateY(rotateAngle[1]);
    transform->RotateZ(rotateAngle[2]);
    transform->Translate(rotateCenter[0], rotateCenter[1], rotateCenter[2]);
    volume->SetUserTransform(transform);
}

/*
 * @Func: this function can transform data type of imagedata;
 * @Paramters: inputImagedata is a input data type of vtkImageData ,
 *             rotateCenter is a double type array ,
 *             rotateAngle is a double type array,
 *             dims is data's dimensions;
 * @Created by cxfeng
 * @return: 0;
 */
void RotationFunction::rotateImageData(vtkImageData *inputImagedata, double *rotateCenter, double *rotateAngle, int *dims){
    vtkImageData* outputImagedata =
            vtkImageData::New ();
    outputImagedata->SetDimensions(inputImagedata->GetDimensions());
    outputImagedata->SetSpacing(inputImagedata->GetSpacing());
    outputImagedata->AllocateScalars(VTK_UNSIGNED_CHAR, 0);

    unsigned char *pTemp = (unsigned char *) outputImagedata->GetScalarPointer();
    unsigned char *ptr = (unsigned char *) inputImagedata->GetScalarPointer();

    vtkSmartPointer<vtkTransform> transform =
            vtkSmartPointer<vtkTransform>::New();
    transform->PostMultiply();
    transform->Translate(-rotateCenter[0], -rotateCenter[1], -rotateCenter[2]);
    transform->RotateX(rotateAngle[0]);
    transform->RotateY(rotateAngle[1]);
    transform->RotateZ(rotateAngle[2]);
    transform->Translate(rotateCenter[0], rotateCenter[1], rotateCenter[2]);

    vtkSmartPointer<vtkMatrix4x4> matrix = transform->GetMatrix();

    int w = 1;
    float inputVector[4]{0};
    float outputVector[4]{0};

    for (int z = 0; z < dims[2]; ++z){
        for (int y = 0; y < dims[1]; ++y){
            for (int x = 0; x < dims[0]; ++x){
                int data = (int)(*(ptr + x + (y*dims[0]) + (z*dims[1]*dims[0])));

                if(data){
                    inputVector[0] = x;
                    inputVector[1] = y;
                    inputVector[2] = z;
                    inputVector[3] = w;
                    matrix->MultiplyPoint(matrix->GetData(), inputVector, outputVector);
                    int outsize = (int)outputVector[0] + ((int)outputVector[1]*dims[0]) + ((int)outputVector[2]*dims[1]*dims[0]);

                    if(outsize >= 0 && outsize < dims[0]*dims[1]*dims[2] ){
                        *(pTemp + (int)outputVector[0] + ((int)outputVector[1]*dims[0]) + ((int)(outputVector[2])*dims[1]*dims[0])) =
                                *(ptr + x + (y*dims[0]) + (z*dims[1]*dims[0]));

                    }
                }
            }
        }
    }
    inputImagedata->ShallowCopy (outputImagedata);
    outputImagedata->Delete ();
}

void RotationFunction::ThereDimensionalReconstruction(){

    readMaskBinary(pelDir, vtkUcharArray_Pel, pelImageData);
    readMaskBinary(leftDir, vtkUcharArray_Left, leftImageData);
    readMaskBinary(rightDir, vtkUcharArray_Right, rightImageData);

#if 0
    clock_t time1 = clock ();
    rotateImageData(pelImageData, pelRotateCenter, pelRotateAngle, dicomDimension);
    time1 = clock () - time1;
    std::cout<<"rotateImageData time :"<<(double)time1 / CLOCKS_PER_SEC<<endl;

    clock_t time2 = clock ();
    rotateImageData(leftImageData, leftRotateCenter, leftRotateAngle, dicomDimension);
    time2 = clock () - time2;
    std::cout<<"rotateImageData time :"<<(double)time2 / CLOCKS_PER_SEC<<endl;

    clock_t time3 = clock ();
    rotateImageData(rightImageData, rightRotateCenter, rightRotateAngle, dicomDimension);
    time3 = clock () - time3;
    std::cout<<"rotateImageData time :"<<(double)time3 / CLOCKS_PER_SEC<<endl;
#endif


    pelVolumeMapper->SetInputData(pelImageData);
    leftVolumeMapper->SetInputData(leftImageData);
    rightVolumeMapper->SetInputData(rightImageData);
    readDicomSeries(dicomDir, dicomImageData);
    dicomVolumeMapper->SetInputData(dicomImageData);

    pelVolume->SetMapper(pelVolumeMapper);
    leftVolume->SetMapper(leftVolumeMapper);
    rightVolume->SetMapper(rightVolumeMapper);
    dicomVolume->SetMapper(dicomVolumeMapper);

#if 1
    clock_t time1 = clock ();
    rotateVolume(pelVolume, pelRotateCenter, pelRotateAngle);
    time1 = clock () - time1;
    std::cout<<"rotateVolume time :"<<(double)time1 / CLOCKS_PER_SEC<<endl;

    clock_t time2 = clock ();
    rotateVolume(leftVolume, leftRotateCenter, leftRotateAngle);
    time2 = clock () - time2;
    std::cout<<"rotateVolume time :"<<(double)time2 / CLOCKS_PER_SEC<<endl;

    clock_t time3 = clock ();
    rotateVolume(rightVolume, rightRotateCenter, rightRotateAngle);
    time3 = clock () - time3;
    std::cout<<"rotateVolume time :"<<(double)time3 / CLOCKS_PER_SEC<<endl;

#endif
    Renderer->AddVolume(pelVolume);
    Renderer->AddVolume(leftVolume);
    Renderer->AddVolume(rightVolume);
//    Renderer->AddVolume(dicomVolume);
    Renderer->SetBackground(0, 0, 0);
    renderWindowInteractor->SetRenderWindow(renderWindow);
    volumeProperty->SetIndependentComponents(independentComponents);
    volumeOpacity->AddSegment(3, 0, 1000, 1);
    color->AddRGBSegment(0.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0);
    volumeProperty->SetColor(color);
    volumeProperty->SetScalarOpacity(volumeOpacity);
    volumeProperty->SetInterpolationTypeToLinear();

    pelVolume->SetProperty(volumeProperty);
    leftVolume->SetProperty(volumeProperty);
    rightVolume->SetProperty(volumeProperty);
    dicomVolume->SetProperty(volumeProperty);

    renderWindow->SetSize(1600, 1000);
    renderWindow->AddRenderer(Renderer);
    renderWindowInteractor->SetRenderWindow(renderWindow);

    coordinateMarkerWidget();
    markRotateCenter(leftRotateCenter);

    renderWindow->Render();
    renderWindowInteractor->Start();
    putPicturesInPdfFile();


}

void RotationFunction::savePngPictures(std::string pngSaveName){
    vtkWindowToImageFilter* WindowToImageFilter =
            vtkWindowToImageFilter::New();
    vtkPNGWriter* PNGWriter =
            vtkPNGWriter::New();
    WindowToImageFilter->SetInput(renderWindow);
    WindowToImageFilter->Update();
    PNGWriter->SetInputConnection(WindowToImageFilter->GetOutputPort());
    PNGWriter->SetFileName(pngSaveName.c_str());
    PNGWriter->SetFilePattern("png");
    PNGWriter->Write();
    PNGWriter->Delete();
    WindowToImageFilter->Delete();
}

void RotationFunction::show(){
    inititializationData();
    ThereDimensionalReconstruction();
//    savePNGPictures();
    return;
}

void RotationFunction::changeCameraFieldOfView(double AzimuthDegree, double ElevationDegree){
    renderWindow->Render();
    Renderer->GetActiveCamera()->Elevation(ElevationDegree);
    renderWindow->Render();
    Renderer->GetActiveCamera()->Azimuth(AzimuthDegree);

}

void RotationFunction::putPicturesInPdfFile(){
    changeCameraFieldOfView(0, 90);
    savePngPictures("1");
    changeCameraFieldOfView(45, 0);
    savePngPictures("2");
    changeCameraFieldOfView(45, 0);
    savePngPictures("3");
    changeCameraFieldOfView(45, 0);
    savePngPictures("4");
    changeCameraFieldOfView(45, 0);
    savePngPictures("5");
    changeCameraFieldOfView(45, 0);
    savePngPictures("6");
    changeCameraFieldOfView(45, 0);
    savePngPictures("7");

    std::cout<<"Here is save PDF file function"<<endl;

}

void RotationFunction::orienntationMarkerWidget(){

        double rgba[4]{0.0, 0.0, 0.0, 0.0};
        widget->SetOutlineColor(rgba[0], rgba[1], rgba[2]);
        widget->SetOrientationMarker( oriAxesActor );
        widget->SetInteractor( renderWindowInteractor );
        widget->SetViewport( 0, 0, 0.3, 0.3 );
        widget->SetEnabled( 1 );
        widget->InteractiveOn();
}

void RotationFunction::coordinateMarkerWidget(){

    oriAxesActor->SetPosition(0, 0, 0);
    oriAxesActor->SetTotalLength(512, 512, 600);
    oriAxesActor->SetShaftType(0);
    oriAxesActor->SetAxisLabels(0);
    oriAxesActor->SetCylinderRadius(0.002);
    Renderer->AddActor(oriAxesActor);
}

void RotationFunction::writeOutBinaryFiles(vtkImageData *inputImagedata, std::string filename){
    char *buffer = new char[dicomSize]{0};
    ofstream myFile(filename.c_str(), ios::out | ios::binary);
    if(!myFile){cout<<"false!!!"<<endl;}
    myFile.write(buffer, dicomSize);
    myFile.close();
}

void RotationFunction::markRotateCenter(double *rotateCenter){
    vtkSmartPointer<vtkSphereSource> pelSphere =
            vtkSmartPointer<vtkSphereSource>::New();
    pelSphere->SetRadius( 10 );
    vtkSmartPointer<vtkSphereSource> leftSphere =
            vtkSmartPointer<vtkSphereSource>::New();
    leftSphere->SetRadius( 10 );
    vtkSmartPointer<vtkSphereSource> rightSphere =
            vtkSmartPointer<vtkSphereSource>::New();
    rightSphere->SetRadius( 10 );

    vtkSmartPointer<vtkPolyDataMapper> pelSphereMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    pelSphereMapper->SetInputConnection( pelSphere->GetOutputPort() );
    vtkSmartPointer<vtkPolyDataMapper> leftSphereMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    leftSphereMapper->SetInputConnection( leftSphere->GetOutputPort() );
    vtkSmartPointer<vtkPolyDataMapper> rightSphereMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    rightSphereMapper->SetInputConnection( rightSphere->GetOutputPort() );

    vtkSmartPointer<vtkActor> pelSphereActor =
            vtkSmartPointer<vtkActor>::New();
    pelSphereActor->SetMapper( pelSphereMapper );
    pelSphereActor->GetProperty ()->SetColor (255, 0 ,0);
    vtkSmartPointer<vtkActor> leftSphereActor =
            vtkSmartPointer<vtkActor>::New();
    leftSphereActor->SetMapper( leftSphereMapper );
    leftSphereActor->GetProperty ()->SetColor (0, 255 ,0);
    vtkSmartPointer<vtkActor> rightSphereActor =
            vtkSmartPointer<vtkActor>::New();
    rightSphereActor->SetMapper( rightSphereMapper );
    rightSphereActor->GetProperty ()->SetColor (0, 0 ,255);

    pelSphereActor->SetPosition(155, 211.188, 83);
    leftSphereActor->SetPosition(340.088, 251.34, 157.613);
    rightSphereActor->SetPosition(163.947, 255.421, 166.661);

    Renderer->AddActor(pelSphereActor);
    Renderer->AddActor(leftSphereActor);
    Renderer->AddActor(rightSphereActor);

}





