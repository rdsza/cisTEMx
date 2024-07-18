#include <cistem_config.h>

#ifdef ENABLEGPU
#include "../../gpu/gpu_core_headers.h"
#include "../../gpu/DeviceManager.h"
#include "../../gpu/TemplateMatchingCore.h"
#else
#include "../../core/core_headers.h"
#endif

#include "../../constants/constants.h"

#if defined(ENABLE_FastFFT) && defined(ENABLEGPU)
#include "../../ext/FastFFT/include/FastFFT.h"
#endif


class
        MatchTemplateApp : public MyApp {
  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );
    void MasterHandleProgramDefinedResult(float* result_array, long array_size, int result_number, int number_of_expected_results);
    void ProgramSpecificInit( );

    // for master collation

    // ArrayOfAggregatedTemplateResults aggregated_results;

    float GetMaxJobWaitTimeInSeconds( ) { return 360.0f; }

  private:
};

IMPLEMENT_APP(MatchTemplateApp)

void MatchTemplateApp::DoInteractiveUserInput( ) {
    wxString input_search_images;
    wxString search_templates;

    wxString cc_output_file; // cross correlation output file 

    float pixel_size              = 1.0f;
    float voltage_kV              = 300.0f;
    float spherical_aberration_mm = 2.7f;
    float amplitude_contrast      = 0.07f;
    float defocus1                = 10000.0f;
    float defocus2                = 10000.0f;

    float defocus_angle;
    float padding = 1.0;

    UserInput* my_input = new UserInput("Perform Cross Correlations on image stack", 1.00);

    // get input
    input_search_images       = my_input->GetFilenameFromUser("Input images to be searched", "", "input.mrc", false);
    search_templates          = my_input->GetFilenameFromUser("Input template stack", "", "input.mrc", false);
    defocus1                  = my_input->GetFloatFromUser("Defocus1 (angstroms)", "Defocus1 for the input image", "10000", 0.0);
    defocus2                  = my_input->GetFloatFromUser("Defocus2 (angstroms)", "Defocus2 for the input image", "10000", 0.0);
    defocus_angle             = my_input->GetFloatFromUser("Defocus Angle (degrees)", "Defocus Angle for the input image", "0.0");
    int first_search_position = 0; //let's try with zero
    int last_search_position  = 91; //maybe 92
    delete my_input;

    my_current_job.ManualSetArguments("ttfffii", input_search_images.ToUTF8( ).data( ), search_templates.ToUTF8( ).data( ),
    defocus1, defocus2, defocus_angle, first_search_position, last_search_position);
  

}

bool MatchTemplateApp::DoCalculation( ) {
    //ends at 1253 in other file
    //Bring inputs over from input function
    wxString input_search_images_filename  = my_current_job.arguments[0].ReturnStringArgument( );
    wxString search_templates_filename     = my_current_job.arguments[1].ReturnStringArgument( );
    float    defocus1                      = my_current_job.arguments[2].ReturnFloatArgument( );
    float    defocus2                      = my_current_job.arguments[3].ReturnFloatArgument( );
    float    defocus_angle                 = my_current_job.arguments[4].ReturnFloatArgument( );
    int      first_search_position         = my_current_job.arguments[5].ReturnIntegerArgument( );
    int      last_search_position          = my_current_job.arguments[6].ReturnIntegerArgument( );
    float pixel_size                       = 1.5f;
    float voltage_kV                       = 300.0f;
    float spherical_aberration_mm          = 2.7f;
    float amplitude_contrast               = 0.07f;
    float padding                          = 1.0;
    //what are these for? 
    // I don't think we need these
    //float    pixel_size_search_range   = 0.1f;
    //float    pixel_size_step           = 0.02f;
    Curve whitening_filter;
    Curve number_of_terms;
    float phase_shift;
    //float defocus_step = 0.0f;
    //float defocus_search_range = 0.0f;
    long pixel_counter;
    int current_search_position;
    int current_x;
    int current_y;
    ImageFile input_search_image_file;
    ImageFile search_templates_file;

    Image input_image;
    Image padded_reference;
    Image search_templates;
    Image template_reconstruction;
    Image current_projection;
    Image padded_projection;
    Image projection_filter;
    Image correlation_pixel_sum_image;
    Image max_intensity_projection;
    double* correlation_pixel_sum            = new double[input_image.real_memory_allocated];
    double* correlation_pixel_sum_of_squares = new double[input_image.real_memory_allocated];
    ZeroDoubleArray(correlation_pixel_sum, input_image.real_memory_allocated);
    ZeroDoubleArray(correlation_pixel_sum_of_squares, input_image.real_memory_allocated);

    
    input_search_image_file.OpenFile(input_search_images_filename.ToStdString( ), false);
    search_templates_file.OpenFile(search_templates_filename.ToStdString( ), false);
    input_image.ReadSlice(&input_search_image_file, 1);
    //Get ctf object to use for calculating ctf
    CTF input_ctf;
    input_ctf.Init(voltage_kV, spherical_aberration_mm, amplitude_contrast, defocus1, defocus2, defocus_angle, 0.0, 0.0, 0.0, pixel_size, deg_2_rad(phase_shift));
    // do Template Match
    //Do we need to do the factorization? 
    //Allocate space
    current_projection.Allocate(search_templates_file.ReturnXSize( ), search_templates_file.ReturnXSize( ), false);
    projection_filter.Allocate(search_templates_file.ReturnXSize( ), search_templates_file.ReturnXSize( ), false);
    template_reconstruction.Allocate(search_templates.logical_x_dimension, search_templates.logical_y_dimension, search_templates.logical_z_dimension, true);
    padded_reference.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    max_intensity_projection.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    max_intensity_projection.SetToConstant(-FLT_MAX);
    padded_reference.SetToConstant(0.0f);
    max_intensity_projection.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    max_intensity_projection.SetToConstant(-FLT_MAX);
    if ( padding != 1.0f )
        padded_projection.Allocate(search_templates_file.ReturnXSize( ) * padding, search_templates_file.ReturnXSize( ) * padding, false);

    // I think we need all this 
    whitening_filter.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
    number_of_terms.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));

    wxDateTime my_time_out;
    wxDateTime my_time_in;
    // remove outliers
    // This won't work for movie frames (13.0 is used in unblur) TODO use poisson stats
    input_image.ReplaceOutliersWithMean(5.0f);
    input_image.ForwardFFT( );
    input_image.SwapRealSpaceQuadrants( );

    input_image.ZeroCentralPixel( );
    input_image.Compute1DPowerSpectrumCurve(&whitening_filter, &number_of_terms);
    whitening_filter.SquareRoot( );
    whitening_filter.Reciprocal( );
    whitening_filter.MultiplyByConstant(1.0f / whitening_filter.ReturnMaximumValue( ));

    input_image.ApplyCurveFilter(&whitening_filter);
    input_image.ZeroCentralPixel( );
    // Note: we are dividing by the sqrt of the sum of squares, so the variance in the images 1/N, not 1. This is where the need to multiply the mips by sqrt(N) comes from.
    // Dividing by sqrt(input_image.ReturnSumOfSquares() / N) would result in a properly normalized CCC value.
    input_image.DivideByConstant(sqrtf(input_image.ReturnSumOfSquares( )));
    //input_image.QuickAndDirtyWriteSlice("/tmp/white.mrc", 1);
    //exit(-1);

    int total_correlation_positions = last_search_position;
    int total_correlation_positions_per_thread = total_correlation_positions;

    ProgressBar* my_progress;

    //Loop over ever search position

    wxPrintf("\n\tFor image id %i\n", image_number_for_gui);
    wxPrintf("Searching %i positions on the Euler sphere (first-last: %i-%i)\n", last_search_position - first_search_position, first_search_position, last_search_position);
    wxPrintf("Searching %i rotations per position.\n", number_of_rotations);
    wxPrintf("There are %li correlation positions total.\n\n", total_correlation_positions);

    wxPrintf("Performing Search...\n\n");

    for ( int current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++ ) {
        // make the projection filter, which will be CTF * whitening filter
        // if defocus step is zero, can we just get rid of the defocus_step and defocus_i?
        input_ctf.SetDefocus((defocus1 + float(defocus_i) * defocus_step) / pixel_size, (defocus2 + float(defocus_i) * defocus_step) / pixel_size, deg_2_rad(defocus_angle));
        //            input_ctf.SetDefocus((defocus1 + 200) / pixel_size, (defocus2 + 200) / pixel_size, deg_2_rad(defocus_angle));
        projection_filter.CalculateCTFImage(input_ctf);
        projection_filter.ApplyCurveFilter(&whitening_filter);
        
        
        if ( padding != 1.0f ) {
            template_reconstruction.ReadSlice(&search_templates_file, current_search_position); //changed to ReadSlice
            padded_projection.SwapRealSpaceQuadrants( );
            padded_projection.BackwardFFT( );
            padded_projection.ClipInto(&current_projection);
            current_projection.ForwardFFT( );
            }
            else {
                template_reconstruction.ReadSlice(&search_templates_file, current_search_position); //  changed to ReadSlice
                current_projection.SwapRealSpaceQuadrants( );
                }
        current_projection.MultiplyPixelWise(projection_filter);
        current_projection.BackwardFFT( );
        //current_projection.ReplaceOutliersWithMean(6.0f);

        current_projection.AddConstant(-current_projection.ReturnAverageOfRealValuesOnEdges( ));

        // We want a variance of 1 in the padded FFT. Scale the small SumOfSquares (which is already divided by n) and then re-divide by N.
        float variance = current_projection.ReturnSumOfSquares( ) * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels - powf(current_projection.ReturnAverageOfRealValues( ) * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels, 2);
        current_projection.DivideByConstant(sqrtf(variance));
        current_projection.ClipIntoLargerRealSpace2D(&padded_reference);

        // Note: The real space variance is set to 1.0 (for the padded size image) and that results in a variance of N in the FFT do to the scaling of the FFT,
        // but the FFT values are divided by 1/N so the variance becomes N / (N^2) = is 1/N
        padded_reference.ForwardFFT( );
        // Zeroing the central pixel is probably not doing anything useful...
        padded_reference.ZeroCentralPixel( );       

#ifdef MKL
                    // Use the MKL
                    vmcMulByConj(padded_reference.real_memory_allocated / 2, reinterpret_cast<MKL_Complex8*>(input_image.complex_values), reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), VML_EP | VML_FTZDAZ_ON | VML_ERRMODE_IGNORE);
#else
                    for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
                          padded_reference.complex_values[pixel_counter] = conj(padded_reference.complex_values[pixel_counter]) * input_image.complex_values[pixel_counter];
                    }
#endif
        // Note: the cross correlation will have variance 1/N (the product of variance of the two FFTs assuming the means are both zero and the distributions independent.)
        // Taking the inverse FFT scales this variance by N resulting in a MIP with variance 1
        padded_reference.BackwardFFT( );
        // update mip, and histogram..
        pixel_counter = 0;
        for ( current_y = 0; current_y < max_intensity_projection.logical_y_dimension; current_y++ ) {
            for ( current_x = 0; current_x < max_intensity_projection.logical_x_dimension; current_x++ ) {
                // first mip
                if ( padded_reference.real_values[pixel_counter] > max_intensity_projection.real_values[pixel_counter] ) {
                    max_intensity_projection.real_values[pixel_counter] = padded_reference.real_values[pixel_counter];
                    }
                    pixel_counter++;
            }
            pixel_counter += padded_reference.padding_jump_value;
            }

            // Write padded_reference to file as slices (92 slices) (imagesize_x,image_size_y,92)

            //                    correlation_pixel_sum.AddImage(&padded_reference);
            for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated; pixel_counter++ ) {
                correlation_pixel_sum[pixel_counter] += padded_reference.real_values[pixel_counter];
                }
                padded_reference.SquareRealValues( );
                //                    correlation_pixel_sum_of_squares.AddImage(&padded_reference);
                for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated; pixel_counter++ ) {
                    correlation_pixel_sum_of_squares[pixel_counter] += padded_reference.real_values[pixel_counter];
                    }

                    //max_intensity_projection.QuickAndDirtyWriteSlice("/tmp/mip.mrc", 1);

                    current_projection.is_in_real_space = false;
                    padded_reference.is_in_real_space   = true;

    }

    //what about code 890 - 1040?
    // write out one single MIP
    //is the mip just a single float value? MIP is an image max_intensity_projection object
}