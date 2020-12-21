/*
This example shows how to perform Hand-Eye calibration.
*/

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

using Zivid.NET.Calibration;
using Duration = Zivid.NET.Duration;

class Program
{
    static void Main()
    {
        try
        {
            var zivid = new Zivid.NET.Application();

            Console.WriteLine("Connecting to camera");
            var camera = zivid.ConnectCamera();
            var inputs = ReadInputs(camera);

            Console.WriteLine("Performing hand-eye calibration");
            var calibrationResult = Calibrator.CalibrateEyeToHand(inputs);

            if (calibrationResult)
            {
                Console.WriteLine("{0}\n{1}\n{2}", "Hand-eye calibration OK", "Result:", calibrationResult);
            }
            else
            {
                Console.WriteLine("Hand-eye calibration FAILED");
                Environment.ExitCode = 1;
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine("Error: {0}", ex.Message);
            Environment.ExitCode = 1;
        }
    }

    static List<HandEyeInput> ReadInputs(Zivid.NET.Camera camera)
    {
        var input = new List<HandEyeInput>();
        var currentPoseId = 0U;
        var beingInput = true;

        Interaction.ExtendInputBuffer(2048);

        do
        {
            switch (Interaction.EnterCommand())
            {
                case CommandType.AddPose:
                    try
                    {
                        var robotPose = Interaction.EnterRobotPose(currentPoseId);
                        using (var frame = Interaction.AssistedCapture(camera))
                        {
                            Console.Write("Detecting checkerboard in point cloud");
                            var result = Detector.DetectFeaturePoints(frame.PointCloud);

                            if (result)
                            {
                                Console.WriteLine("OK");
                                input.Add(new HandEyeInput(robotPose, result));
                                ++currentPoseId;
                            }
                            else
                            {
                                Console.WriteLine("FAILED");
                            }
                        }
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("Error: {0}", ex.Message);
                        continue;
                    }
                    break;

                case CommandType.Calibrate: beingInput = false; break;

                case CommandType.Unknown: Console.WriteLine("Error: Unknown command"); break;
            }
        } while (beingInput);
        return input;
    }
}

enum CommandType
{
    AddPose,
    Calibrate,
    Unknown
}

class Interaction
{
    // Console.ReadLine only supports reading 256 characters, by default. This limit is modified
    // by calling ExtendInputBuffer with the maximum length of characters to be read.
    public static void ExtendInputBuffer(int size)
    {
        Console.SetIn(new StreamReader(Console.OpenStandardInput(), Console.InputEncoding, false, size));
    }

    public static CommandType EnterCommand()
    {
        Console.Write("Enter command, p (to add robot pose) or c (to perform calibration): ");
        var command = Console.ReadLine().ToLower();

        switch (command)
        {
            case "p": return CommandType.AddPose;
            case "c": return CommandType.Calibrate;
            default: return CommandType.Unknown;
        }
    }

    public static Pose EnterRobotPose(ulong index)
    {
        var elementCount = 16;
        Console.WriteLine(
            "Enter pose with id (a line with {0} space separated values describing 4x4 row-major matrix) : {1}",
            elementCount,
            index);
        var input = Console.ReadLine();

        var elements = input.Split().Where(x => !string.IsNullOrEmpty(x.Trim())).Select(x => float.Parse(x)).ToArray();

        var robotPose = new Pose(elements);
        Console.WriteLine("The following pose was entered: \n{0}", robotPose);
        return robotPose;
    }

    public static Zivid.NET.Frame AssistedCapture(Zivid.NET.Camera camera)
    {
        var suggestSettingsParameters = new Zivid.NET.CaptureAssistant.SuggestSettingsParameters
        {
            AmbientLightFrequency =
                Zivid.NET.CaptureAssistant.SuggestSettingsParameters.AmbientLightFrequencyOption.none,
            MaxCaptureTime = Duration.FromMilliseconds(800)
        };
        var settings = Zivid.NET.CaptureAssistant.Assistant.SuggestSettings(camera, suggestSettingsParameters);
        return camera.Capture(settings);
    }
}
