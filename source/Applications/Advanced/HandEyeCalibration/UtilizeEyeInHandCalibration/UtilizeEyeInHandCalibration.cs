/*
This example shows how to utilize the result of Eye-in-Hand calibration to transform either a (picking) point
coordinates or the entire point cloud from the camera frame to the robot base frame. The YAML files for this sample can
be found under the main instructions for Zivid samples.
*/

using System;
using System.IO;

using YamlDotNet.RepresentationModel;
using MathNet.Numerics.LinearAlgebra;

class Program
{
    static void Main()
    {
        try
        {
            var zivid = new Zivid.NET.Application();

            Console.WriteLine("Reading camera pose in end-effector frame (result of eye-in-hand calibration)");
            var eyeInHandTransformationFile = Environment.GetFolderPath(Environment.SpecialFolder.CommonApplicationData)
                               + "/Zivid/EyeInHandTransform.yaml";

            Console.WriteLine("Reading end-effector pose in robot base frame");
            var endEffectorPoseFile = Environment.GetFolderPath(Environment.SpecialFolder.CommonApplicationData)
                               + "/Zivid/RobotTransform.yaml";

            var eyeInHandTransformation = readTransform(eyeInHandTransformationFile);
            var endEffectorPose = readTransform(endEffectorPoseFile);

            // Converting to Mat.NET matrices for easier computation
            Matrix<float> eyeInHandMatrix = zividToMathDotNet(eyeInHandTransformation);
            Matrix<float> endEffectorMatrix = zividToMathDotNet(endEffectorPose);

            var dataFile = Environment.GetFolderPath(Environment.SpecialFolder.CommonApplicationData) + "/Zivid/ZividGem.zdf";
            Console.WriteLine("Reading ZDF frame from file: " + dataFile);
            var frame = new Zivid.NET.Frame(dataFile);
            var pointCloud = frame.PointCloud;

            Console.WriteLine("Computing camera pose in robot base frame");
            Matrix<float> transformBaseToCamera = endEffectorMatrix * eyeInHandMatrix;

            switch (Interaction.EnterCommand())
            {
                case CommandType.TransformSinglePoint:

                    Console.WriteLine("Transforming single point");

                    // The (picking) point is defined as image coordinates in camera frame. It is hard-coded for the
                    // ZividGem.zdf
                    const int imageCoordinateX = 1357;
                    const int imageCoordinateY = 666;
                    var xyz = pointCloud.CopyPointsXYZW();

                    var pickingPoint = new float[]{xyz[imageCoordinateY, imageCoordinateX, 0],
                                                            xyz[imageCoordinateY, imageCoordinateX, 1],
                                                            xyz[imageCoordinateY, imageCoordinateX, 2],
                                                            xyz[imageCoordinateY, imageCoordinateX, 3]};

                    Vector<float> pointInCameraFrame = CreateVector.DenseOfArray(pickingPoint);
                    Console.WriteLine("Point coordinates in camera frame: " + pointInCameraFrame);

                    Console.WriteLine("Transforming (picking) point from camera to robot base frame");
                    Vector<float> pointInBaseFrame = transformBaseToCamera * pointInCameraFrame;

                    Console.WriteLine("Point coordinates in robot base frame: " + pointInBaseFrame);
                    break;

                case CommandType.TransformPointCloud:

                    Console.WriteLine("Transforming point cloud");
                    var transformMatrix = mathDotNetToZivid(transformBaseToCamera);
                    pointCloud.Transform(transformMatrix);

                    var dataFileTransformed = "ZividGemTransformed.zdf";
                    Console.WriteLine("Saving frame to file: " + dataFileTransformed);
                    frame.Save(dataFileTransformed);
                    break;

                case CommandType.Unknown:
                    Console.WriteLine("Error: Unknown command");
                    break;
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine("Error: {0}", ex.Message);
            Environment.ExitCode = 1;
        }
    }
    static float[,] readTransform(string transformFile)
    {
        using (var reader = new StreamReader(transformFile))
        {
            var yaml = new YamlStream();
            yaml.Load(reader);

            var poseStateNode = (YamlMappingNode)yaml.Documents[0].RootNode;
            var poseStateData = poseStateNode["PoseState"]["data"].ToString();
            var matrixAsString = poseStateData.Trim('[', ']').Split(',');

            int rows = int.Parse(poseStateNode["PoseState"]["rows"].ToString());
            int cols = int.Parse(poseStateNode["PoseState"]["cols"].ToString());

            float[,] zividMatrix = new float[rows, cols];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    zividMatrix[i, j] = float.Parse(matrixAsString[i * rows + j]);
                }
            }
            return zividMatrix;
        }
    }
    static Matrix<float> zividToMathDotNet(float[,] zividMatrix)
    {
        var mathNetMatrix = CreateMatrix.DenseOfArray(zividMatrix);
        return mathNetMatrix;
    }

    static float[,] mathDotNetToZivid(Matrix<float> mathNetMatrix)
    {
        float[,] zividMatrix = mathNetMatrix.ToArray();
        return zividMatrix;
    }
}

enum CommandType
{
    TransformSinglePoint,
    TransformPointCloud,
    Unknown
}

class Interaction
{
    public static CommandType EnterCommand()
    {
        Console.Write("Enter command, s (to transform single point) or p (to transform point cloud): ");
        var command = Console.ReadLine().ToLower();

        switch (command)
        {
            case "s": return CommandType.TransformSinglePoint;
            case "p": return CommandType.TransformPointCloud;
            default: return CommandType.Unknown;
        }
    }
}
