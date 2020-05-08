using System;
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

            Console.WriteLine("Configuring 2D settings");
            // Note: The Zivid SDK supports 2D captures with a single acquisition only
            var settings2D = new Zivid.NET.Settings2D
            {
                Acquisitions = { new Zivid.NET.Settings2D.Acquisition{
                    Aperture = 2.83, ExposureTime = Duration.FromMicroseconds(10000), Gain = 1.0, Brightness = 1.0 } },
                Processing = { Color = { Balance = { Red = 1.0, Blue = 1.0, Green = 1.0 } } }
            };

            Console.WriteLine("Capturing 2D frame");
            using (var frame2D = camera.Capture(settings2D))
            {
                Console.WriteLine("Getting RGBA image");
                var image = frame2D.ImageRGBA();

                var pixelRow = 100;
                var pixelCol = 50;

                var pixelArray = image.ToArray();
                Console.WriteLine("Height: {0}, Width: {1}", pixelArray.GetLength(0), pixelArray.GetLength(1));
                Console.WriteLine("Color at pixel ({0},{1}):  R:{2}  G:{3}  B:{4}  A:{5}",
                                  pixelRow,
                                  pixelCol,
                                  pixelArray[pixelRow, pixelCol].r,
                                  pixelArray[pixelRow, pixelCol].g,
                                  pixelArray[pixelRow, pixelCol].b,
                                  pixelArray[pixelRow, pixelCol].a);

                // Get 3D array of bytes
                var nativeArray = image.ToByteArray();
                Console.WriteLine("Image Height: {0}, Image Width: {1}, Channels: {2}",
                                  nativeArray.GetLength(0),
                                  nativeArray.GetLength(1),
                                  nativeArray.GetLength(2));
                Console.WriteLine("Color at pixel ({0},{1}):  R:{2}  G:{3}  B:{4}  A:{5}",
                                  pixelRow,
                                  pixelCol,
                                  nativeArray[pixelRow, pixelCol, 0],
                                  nativeArray[pixelRow, pixelCol, 1],
                                  nativeArray[pixelRow, pixelCol, 2],
                                  nativeArray[pixelRow, pixelCol, 3]);

                var resultFile = "image.png";
                Console.WriteLine("Saving the image to {0}", resultFile);
                image.Save(resultFile);
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine("Error: " + ex.Message);
            Environment.ExitCode = 1;
        }
    }
}
