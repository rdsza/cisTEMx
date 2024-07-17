#include "../../core/core_headers.h"

class
        fft_write : public MyApp {

            public:
            bool DoCalculation( );
            void DoInteractiveUserInput( );

            private:
        };

IMPLEMENT_APP(fft_write)

void fft_write::DoInteractiveUserInput( ){

    int new_z_size = 1;

    UserInput* my_input = new UserInput("fft_write", 1.0);

    std::string input_filename = my_input->GetFilenameFromUser("Input image file name", "Filename of input image", "input.mrc", true);
    std::string output_filename    = my_input->GetFilenameFromUser("Output image file name", "Filename of output image", "output.mrc", false);

    delete my_input;

    my_current_job.Reset(6);
    my_current_job.ManualSetArguments("tt", input_filename.c_str( ), output_filename.c_str( ));
}

bool fft::DoCalculation( ) {

    std::string input_filename     = my_current_job.arguments[0].ReturnStringArgument( );
    std::string output_filename    = my_current_job.arguments[1].ReturnStringArgument( );
    float       pixel_size;

    ImageFile my_input_file(input_filename, false);
    MRCFile   my_output_file(output_filename, true);

    Image                 my_image;
    Image                 my_amplitude_spectrum;
    // pixel size could be non-square/cubic but we will ignore this here and assume it is square/cubic
    pixel_size = my_input_file.ReturnPixelSize( );

    wxPrintf("\nFourier transforming Images...\n\n");
    ProgressBar* my_progress = new ProgressBar(my_input_file.ReturnNumberOfSlices( ));

    for ( long image_counter = 0; image_counter < my_input_file.ReturnNumberOfSlices( ); image_counter++ ) {
        my_image.ReadSlice(&my_input_file, image_counter + 1);

        my_amplitude_spectrum.Allocate(my_image.logical_x_dimension, my_image.logical_y_dimension, true);

        my_image.ForwardFFT( );
        my_image.ComputeAmplitudeSpectrumFull2D(&my_amplitude_spectrum);
        my_amplitude_spectrum.UpdateDistributionOfRealValues(&my_distribution);

        my_amplitude_spectrum.WriteSlice(&my_output_file, image_counter + 1);
        // next line is buggy because only the last image will determine image header stats
        my_progress->Update(image_counter + 1);
    }

    delete my_progress;
    wxPrintf("\n\n");

    float std = my_distribution.GetSampleVariance( );
    if ( std > 0.0 ) {
        std = sqrt(std);
    }
    
    my_output_file.SetDensityStatistics(my_distribution.GetMinimum( ), my_distribution.GetMaximum( ), my_distribution.GetSampleMean( ), std);

    my_output_file.SetPixelSize(pixel_size);
    my_output_file.WriteHeader( );

    return true;
}
