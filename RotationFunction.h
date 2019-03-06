//
// Created by cxfeng on 19-2-21.
//

#ifndef VOLUMERENDERING_ROTATIONFUNCTION_H
#define VOLUMERENDERING_ROTATIONFUNCTION_H
/*
 * Qt Matrix
 */

/*
 * VTK read or source includes
 */
#include "vtkDICOMImageReader.h"
#include <vtkSphereSource.h>

/*
 * VTK data structure includes
 */
#include <vtkUnsignedCharArray.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>

/*
 * VTK function includes
 */
#include <vtkSmartPointer.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkMatrix4x4.h>
#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

/*
 * VTK Property
 */
#include <vtkProperty.h>
#include "vtkVolumeProperty.h"
#include <vtkTextProperty.h>

/*
 * VTK rendering engine includes
 */
#include <vtkPolyDataMapper.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include "vtkVolume.h"
#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkTextActor.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

/*
 * C++ lib
 */
#include <time.h>
#include <fstream>
#include <unistd.h>
/*
 * Class RotationFunction aims is transform imagedata and volume;
 * Created by cxfeng;
 */
class RotationFunction {

public:
    explicit RotationFunction();

    ~RotationFunction();

    void show();

    struct RotateParam
    {
        double rotateAngle[3];
        double rotateCenter[3];
    };

private:
    void readMaskBinary(std::string dataDir, vtkUnsignedCharArray *ptr1, vtkImageData *ptr);

    void readDicomSeries(std::string dataDir, vtkImageData *ptr);

    void maskAndDicomAreMultipliedToDivide(vtkImageData *imageData1, vtkImageData *imageData2);

    void maskAndDicomAreMultipliedToDivideWhenYAxisMirror(vtkImageData *imageData1, vtkImageData *imageData2);

    void rotateVolume(vtkVolume *volume, double rotateCenter[3], double rotateAngle[3]);

    void rotateImageData(vtkImageData *inputImagedata, double *rotateCenter, double *rotateAngle, int *dims);

    void writeOutBinaryFiles(vtkImageData *inputImagedata, std::string filename);

    void ThereDimensionalReconstruction();

    void changeCameraFieldOfView(double AzimuthDegree, double ElevationDegree);

    void savePngPictures(std::string pngSaveName);

    void putPicturesInPdfFile();

    void orienntationMarkerWidget();

    void coordinateMarkerWidget();

    void markRotateCenter(double *rotateCenter);

    std::string pelDir;
    std::string leftDir;
    std::string rightDir;
    std::string dicomDir;
    int *dicomDimension;
    double *spacings;
    const int dicomSize = 512*512*600;
    std::string pngSaveName;

    double *pelRotateAngle;
    double *pelRotateCenter;
    double *leftRotateAngle;
    double *leftRotateCenter;
    double *rightRotateAngle;
    double *rightRotateCenter;

    double *cameraPosition;
    double *cameraViewUp;
    double *cameraFocalPoint;
    double *cameraClippingRange;

    bool independentComponents;

    vtkSmartPointer<vtkUnsignedCharArray> vtkUcharArray_Pel;
    vtkSmartPointer<vtkUnsignedCharArray> vtkUcharArray_Left;
    vtkSmartPointer<vtkUnsignedCharArray> vtkUcharArray_Right;

    // 4 imagedata
    vtkSmartPointer<vtkImageData> pelImageData;
    vtkSmartPointer<vtkImageData> leftImageData;
    vtkSmartPointer<vtkImageData> rightImageData;
    vtkSmartPointer<vtkImageData> dicomImageData;

    // 4 mapper
    vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> pelVolumeMapper;
    vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> leftVolumeMapper;
    vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> rightVolumeMapper;
    vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> dicomVolumeMapper;

    // 4 volume
    vtkSmartPointer<vtkVolume> pelVolume;
    vtkSmartPointer<vtkVolume> leftVolume;
    vtkSmartPointer<vtkVolume> rightVolume;
    vtkSmartPointer<vtkVolume> dicomVolume;

    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> interactorStyle;
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
    vtkSmartPointer<vtkVolumeProperty> volumeProperty;
    vtkSmartPointer<vtkPiecewiseFunction> volumeOpacity;
    vtkSmartPointer<vtkColorTransferFunction> color;
    vtkSmartPointer<vtkRenderer> Renderer;

    vtkSmartPointer<vtkAxesActor> oriAxesActor;
    vtkSmartPointer<vtkOrientationMarkerWidget> widget;
    vtkSmartPointer<vtkCamera> camera;

    void inititializationData();

};


#endif //VOLUMERENDERING_ROTATIONFUNCTION_H
