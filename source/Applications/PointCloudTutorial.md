
## Introduction

This tutorial describes how to use Zivid SDK to work with [Point Cloud][kb-point_cloud-url] data.

1. [Frame](#frame)
   1. [Capture](#capture)
   2. [Load](#load)
   3. [Visualize](#visualize)
3.  [Point Cloud](#point-cloud)
    1. [Get handle](#get-handle)
    2. [Copy](#copy)
    3.  [Transform](#transform)
    4. [Downsample](#downsample)

### Prerequisites

You should have installed Zivid SDK and cloned C# samples. For more details see [Instructions][installation-instructions-url].

## Frame

### Capture

The ```Zivid.NET.Frame``` can be captured with a camera ([go to source][frame-capture]).
```csharp
var frame = camera.Capture(settings);
```
Check [Capture Tutorial][capture-tutorial] for more details.

### Load
It can also be loaded from a ZDF file ([go to source][frame-from-file]).
```csharp
var frame = new Zivid.NET.Frame("/path/to/Zivid3D.zdf");
```
### Visualize

Having the frame allows you to visualize the point cloud ([go to source][visualize-point-cloud]).

```csharp
var visualizer = new Zivid.NET.Visualization.Visualizer();

visualizer.Show(frame);
visualizer.ShowMaximized();
visualizer.ResetToFit();

visualizer.Run();
```

## Point Cloud

### Get handle

You can now get a handle to the point cloud data on the GPU ([go to source][point-cloud]).
```csharp
var pointCloud = frame.PointCloud;
```
Getting the property ```Zivid.NET.Frame.PointCloud``` does not perform any copying from GPU memory.

Note: ```Zivid.NET.Camera.Capture()``` method returns as soon as the camera is done capturing. At that point the handle from ```Zivid.NET.Frame.PointCloud``` is available instantly as well. However, the actual point cloud data becomes available only after the GPU is finished processing. This processing will occur on the GPU in the background, and any calls to data-copy methods (see section [Copy](#copy) below) will wait for this to finish before proceeding with the requested copy operation.

### Copy

You can now selectively copy data based on what is required. This is the complete list of output data formats and how to copy them from the GPU.


|Return type|Methods for copying from GPU|Data per pixel|Total data copied|
|-|-|-|-|
|```float[height,width,3]```| ```PointCloud.CopyPointsXYZ()```| 12 bytes |28 MB |
|```float[height,width,4]```| ```PointCloud.CopyPointsXYZW()```| 16 bytes |37 MB |
|```float[height,width,1]```| ```PointCloud.CopyPointsZ()```| 4 bytes |9 MB |
|```byte[height,width,4]```| ```PointCloud.CopyColorsRGBA()```| 4 bytes |9 MB |
|```Zivid.NET.ImageRGBA```| ```PointCloud.CopyImageRGBA()```| 4 bytes |9 MB |
|```Zivid.NET.PointXYZColorRGBA[height, width]```| ```PointCloud.CopyPointsXYZColorsRGBA()```| 16 bytes |37 MB |
|```Zivid.NET.PointXYZColorBGRA[height, width]```| ```PointCloud.CopyPointsXYZColorsBGRA()```| 16 bytes |37 MB |
|```float[height,width]```| ```PointCloud.CopySNRs()```| 4 bytes |9 MB |

#### Copy selected data from GPU to system memory (Zivid-allocated)

If you are only concerned about e.g. RGB color data of the point cloud, you can copy only that data to the system memory ([go to source][copy]).
```csharp
var colorsRGBA = frame.PointCloud.CopyColorsRGBA(); /* Colors are copied from the GPU and into a
three-dimensional array of bytes (byte[height,width,4]), representing a 4-channel RGBA image. */
```

## Transform

You may want to change the point cloud's origin from the camera to the robot base frame or scale the point cloud to e.g. change it from millimeters to meters. This can be done by transforming the point cloud using ```Zivid.NET.PointCloud.Transform(transformationMatrix)```.

```csharp
pointCloud.Transform(transformationMatrix);
```

## Downsample


Sometimes you might not need as dense point cloud as given from the camera and you want to downsample the point cloud. Downsampling can be done in-place by using ```Zivid.NET.PointCloud.Downsample(Zivid.NET.PointCloud.Downsampling)``` ([go to source][downsample]).

```csharp
pointCloud.Downsample(Zivid.NET.PointCloud.Downsampling.By2x2);
```
Another way to downsample a point cloud is by using ```Zivid.NET.PointCloud.Downsampled(Zivid.NET.PointCloud.Downsampling)```. This does not modify the current point cloud but returns the downsampled point cloud as a new point cloud instance. 

```csharp
downsampledPointCloud = pointCloud.Downsampled(Zivid.NET.PointCloud.Downsampling.By2x2);
```

Supported downsampling rates are: by2x2, by3x3 and by4x4.

## Conclusion

This tutorial shows how to use the Zivid SDK to extract the point cloud, manipulate it, transform it, and visualize it.

[//]: ### "Recommended further reading"

[installation-instructions-url]: ../../README.md#instructions
[frame-from-file]:Basic/FileFormats/ReadIterateZDF/ReadIterateZDF.cs#L15-L18
[frame-capture]:../Camera/Basic/Capture/Capture.cs#L28
[capture-tutorial]:../Camera/Basic/CaptureTutorial.md#L158
[point-cloud]:Basic/FileFormats/ReadIterateZDF/ReadIterateZDF.cs#L21
[copy]:Basic/FileFormats/ReadIterateZDF/ReadIterateZDF.cs#L22-L23
[visualize-point-cloud]:Basic/Visualization/CaptureVis3D/CaptureVis3D.cs#L25-L34
[kb-point_cloud-url]: https://support.zivid.com/latest/rst/reference-articles/zivid-3d-camera-technology/point-cloud-structure-and-output-formats.html
[downsample]:Advanced/Downsample/Downsample.cs#L26]