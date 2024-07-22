#include <cistem_config.h>
#include <iostream>
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
    
    void DoInteractiveUserInput( );
    bool DoCalculation( );
    
    // void MasterHandleProgramDefinedResult(float* result_array, long array_size, int result_number, int number_of_expected_results);
    // void ProgramSpecificInit( );

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
    // input_search_images       = my_input->GetFilenameFromUser("Input images to be searched", "", "input.mrc", false);
    // search_templates          = my_input->GetFilenameFromUser("Input template stack", "", "input.mrc", false);
    // defocus1                  = my_input->GetFloatFromUser("Defocus1 (angstroms)", "Defocus1 for the input image", "10000", 0.0);
    // defocus2                  = my_input->GetFloatFromUser("Defocus2 (angstroms)", "Defocus2 for the input image", "10000", 0.0);
    // defocus_angle             = my_input->GetFloatFromUser("Defocus Angle (degrees)", "Defocus Angle for the input image", "0.0");
    input_search_images       = "/home/useradmin/Match_PCA_template_repo/cisTEM/src/programs/match_pca_template/00040_3_0.mrc";
    search_templates          = "/home/useradmin/Match_PCA_template_repo/cisTEM/src/programs/match_pca_template/1.5_psi_92_peaks.mrc";
    defocus1                  = 13850;
    defocus2                  = 13272;
    defocus_angle             = -4.5;
    cc_output_file        = my_input->GetFilenameFromUser("Output cross correlation filename","", "output_cc.mrc", false);
    // int first_search_position = 0; //let's try with zero
    // int last_search_position  = 91; //maybe 92
    delete my_input;

    // my_current_job.ManualSetArguments("ttfffiit", input_search_images.ToUTF8( ).data( ), search_templates.ToUTF8( ).data( ),
    // defocus1, defocus2, defocus_angle, first_search_position, last_search_position, cc_output_file.ToUTF8( ).data( ));
    my_current_job.ManualSetArguments("t", cc_output_file.ToUTF8( ).data( ));
}

bool MatchTemplateApp::DoCalculation( ) {
    //ends at 1253 in other file
    //Bring inputs over from input function
    wxPrintf("can I print?");
    wxDateTime start_time = wxDateTime::Now( );

    // wxString input_search_images_filename  = my_current_job.arguments[0].ReturnStringArgument( );
    // wxString search_templates_filename     = my_current_job.arguments[1].ReturnStringArgument( );
    // float    defocus1                      = my_current_job.arguments[2].ReturnFloatArgument( );
    // float    defocus2                      = my_current_job.arguments[3].ReturnFloatArgument( );
    // float    defocus_angle                 = my_current_job.arguments[4].ReturnFloatArgument( );
    // int      first_search_position         = my_current_job.arguments[5].ReturnIntegerArgument( );
    // int      last_search_position          = my_current_job.arguments[6].ReturnIntegerArgument( );
    wxString cc_output_filename  = my_current_job.arguments[0].ReturnStringArgument( );
    wxString input_search_images_filename       = "/home/useradmin/Match_PCA_template_repo/cisTEM/src/programs/match_pca_template/00040_3_0.mrc";
    wxString search_templates_filename          = "/home/useradmin/Match_PCA_template_repo/cisTEM/src/programs/match_pca_template/1.5_psi_92_peaks.mrc";
    float defocus1                  = 13850;
    float defocus2                  = 13272;
    float defocus_angle             = -4.5;
   

    int first_search_position = 1; //let's try with zero
    int last_search_position  = 92; //maybe 92
    float pixel_size                       = 1.5f;
    float voltage_kV                       = 300.0f;
    float spherical_aberration_mm          = 2.7f;
    float amplitude_contrast               = 0.07f;
    float padding                          = 1.0f;
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
    // Image max_intensity_projection;
    Image cc_output;
    // double* correlation_pixel_sum            = new double[input_image.real_memory_allocated];
    // double* correlation_pixel_sum_of_squares = new double[input_image.real_memory_allocated];
    // ZeroDoubleArray(correlation_pixel_sum, input_image.real_memory_allocated);
    // ZeroDoubleArray(correlation_pixel_sum_of_squares, input_image.real_memory_allocated);
    
    
    input_search_image_file.OpenFile(input_search_images_filename.ToStdString( ), false);
    search_templates_file.OpenFile(search_templates_filename.ToStdString( ), false);
    input_image.ReadSlice(&input_search_image_file, 1);




 // Resize input image to be factorizable by small numbers
    int factorizable_x;
    int factorizable_y;
    int factor_result_pos;
    int factor_result_neg;
    long   original_input_image_x;
    long   original_input_image_y;
    original_input_image_x = input_image.logical_x_dimension;
    original_input_image_y = input_image.logical_y_dimension;
    factorizable_x         = input_image.logical_x_dimension;
    factorizable_y         = input_image.logical_y_dimension;
    int i;
    float  factor_score;

    bool      DO_FACTORIZATION                       = true;
    bool      MUST_BE_POWER_OF_TWO                   = false; // Required for half-precision xforms
    bool      MUST_BE_FACTOR_OF_FOUR                 = true; // May be faster
    const int max_number_primes                      = 6;
    int       primes[max_number_primes]              = {2, 3, 5, 7, 9, 13};
    float     max_reduction_by_fraction_of_reference = 0.000001f; // FIXME the cpu version is crashing when the image is reduced, but not the GPU
    float     max_increas_by_fraction_of_image       = 0.1f;
    int       max_padding                            = 0; // To restrict histogram calculation
    float     histogram_padding_trim_rescale; // scale the counts to

    // for 5760 this will return
    // 5832 2     2     2     3     3     3     3     3     3 - this is ~ 10% faster than the previous solution BUT
    if ( DO_FACTORIZATION ) {
        for ( i = 0; i < max_number_primes; i++ ) {

            factor_result_neg = ReturnClosestFactorizedLower(original_input_image_x, primes[i], true, MUST_BE_FACTOR_OF_FOUR);
            factor_result_pos = ReturnClosestFactorizedUpper(original_input_image_x, primes[i], true, MUST_BE_FACTOR_OF_FOUR);

            if ( (float)(original_input_image_x - factor_result_neg) < (float)search_templates_file.ReturnXSize( ) * max_reduction_by_fraction_of_reference ) {
                factorizable_x = factor_result_neg;
                break;
            }
            if ( (float)(-original_input_image_x + factor_result_pos) < (float)input_image.logical_x_dimension * max_increas_by_fraction_of_image ) {
                factorizable_x = factor_result_pos;
                break;
            }
        }
        factor_score = FLT_MAX;
        for ( i = 0; i < max_number_primes; i++ ) {

            factor_result_neg = ReturnClosestFactorizedLower(original_input_image_y, primes[i], true, MUST_BE_FACTOR_OF_FOUR);
            factor_result_pos = ReturnClosestFactorizedUpper(original_input_image_y, primes[i], true, MUST_BE_FACTOR_OF_FOUR);

            if ( (float)(original_input_image_y - factor_result_neg) < (float)search_templates_file.ReturnYSize( ) * max_reduction_by_fraction_of_reference ) {
                factorizable_y = factor_result_neg;
                break;
            }
            if ( (float)(-original_input_image_y + factor_result_pos) < (float)input_image.logical_y_dimension * max_increas_by_fraction_of_image ) {
                factorizable_y = factor_result_pos;
                break;
            }
        }
        if ( factorizable_x - original_input_image_x > max_padding )
            max_padding = factorizable_x - original_input_image_x;
        if ( factorizable_y - original_input_image_y > max_padding )
            max_padding = factorizable_y - original_input_image_y;

        if ( ReturnThreadNumberOfCurrentThread( ) == 0 ) {
            wxPrintf("old x, y = %i %i\n  new x, y = %i %i\n", input_image.logical_x_dimension, input_image.logical_y_dimension, factorizable_x, factorizable_y);
        }

        input_image.Resize(factorizable_x, factorizable_y, 1, input_image.ReturnAverageOfRealValuesOnEdges( ));

        search_templates.ReadSlices(&search_templates_file, 1, search_templates_file.ReturnNumberOfSlices( ));
        if ( padding != 1.0f ) {
            search_templates.Resize(search_templates.logical_x_dimension * padding, search_templates.logical_y_dimension * padding, search_templates.logical_z_dimension * padding, search_templates.ReturnAverageOfRealValuesOnEdges( ));
        }


    }


    //Get ctf object to use for calculating ctf
    CTF input_ctf;
    input_ctf.Init(voltage_kV, spherical_aberration_mm, amplitude_contrast, defocus1, defocus2, defocus_angle, 0.0, 0.0, 0.0, pixel_size, deg_2_rad(phase_shift));
    // do Template Match
    //Do we need to do the factorization? 
    //Allocate space
    
    current_projection.Allocate(search_templates_file.ReturnXSize( ), search_templates_file.ReturnXSize( ), false);
    projection_filter.Allocate(search_templates_file.ReturnXSize( ), search_templates_file.ReturnXSize( ), false);
    template_reconstruction.Allocate(search_templates.logical_x_dimension, search_templates.logical_y_dimension, 1, true);
    padded_reference.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    padded_reference.SetToConstant(0.0f);
    // max_intensity_projection.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    // max_intensity_projection.SetToConstant(-FLT_MAX);
    if ( padding != 1.0f )
        padded_projection.Allocate(search_templates_file.ReturnXSize( ) * padding, search_templates_file.ReturnXSize( ) * padding, false);
    cc_output.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);

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
    //may have to change some of these because some of these variables are not used
    // wxPrintf("\n\tFor image id %i\n", image_number_for_gui);
    // wxPrintf("\n\tFor image id %i\n", "00040_3.mrc");
    // wxPrintf("Searching %i positions on the Euler sphere (first-last: %i-%i)\n", last_search_position - first_search_position, first_search_position, last_search_position);
    // wxPrintf("Searching %i rotations per position.\n", number_of_rotations);
    wxPrintf("There are %li correlation positions total.\n\n", total_correlation_positions);

    wxPrintf("Performing Search...\n\n");
        
               //    template_reconstruction.ForwardFFT();
    

        input_ctf.SetDefocus(defocus1, defocus2, deg_2_rad(defocus_angle));
        //            input_ctf.SetDefocus((defocus1 + 200) / pixel_size, (defocus2 + 200) / pixel_size, deg_2_rad(defocus_angle));
        projection_filter.CalculateCTFImage(input_ctf);
        projection_filter.ApplyCurveFilter(&whitening_filter);
      
        // search_templates.ChangePixelSize(&template_reconstruction, 1.5, 0.001f, true);
        template_reconstruction.ZeroCentralPixel( );
        template_reconstruction.SwapRealSpaceQuadrants( );
  wxPrintf("There are %i xs, %i ys, %i zs \n", template_reconstruction.logical_x_dimension, template_reconstruction.logical_y_dimension, template_reconstruction.logical_z_dimension);
        


    for ( int current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++ ) {
        // make the projection filter, which will be CTF * whitening filter
        // if defocus step is zero, can we just get rid of the defocus_step and defocus_i?
        // input_ctf.SetDefocus((defocus1 + float(defocus_i) * defocus_step) / pixel_size, (defocus2 + float(defocus_i) * defocus_step) / pixel_size, deg_2_rad(defocus_angle));
        // wxPrintf("here ");
        // wxPrintf(&search_templates_file);
          //wxPrintf("There are %i xs, %i ys, %i zs \n", template_reconstruction.logical_x_dimension, template_reconstruction.logical_y_dimension, template_reconstruction.logical_z_dimension);

        if ( padding != 1.0f ) {
            
            template_reconstruction.ReadSlice(&search_templates_file, current_search_position); //changed to ReadSlice
            padded_projection.SwapRealSpaceQuadrants( );
            padded_projection.BackwardFFT( );
            padded_projection.ClipInto(&current_projection);
            current_projection.ForwardFFT( );
            }
            else {
                 // wxPrintf("There are %i xs, %i ys, %i zs \n", template_reconstruction.logical_x_dimension, template_reconstruction.logical_y_dimension, template_reconstruction.logical_z_dimension);
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
        for ( current_y = 0; current_y < cc_output.logical_y_dimension; current_y++ ) {
            for ( current_x = 0; current_x < cc_output.logical_x_dimension; current_x++ ) {
                // first mip

                    cc_output.real_values[pixel_counter] = padded_reference.real_values[pixel_counter];
                    
                    pixel_counter++;
            }

            pixel_counter += padded_reference.padding_jump_value;
            }
            cc_output.QuickAndDirtyWriteSlice(cc_output_filename.ToStdString( ), current_search_position + 1);

            // Write padded_reference to file as slices (92 slices) (imagesize_x,image_size_y,92)

    //         //                    correlation_pixel_sum.AddImage(&padded_reference);
    //         for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated; pixel_counter++ ) {
    //             correlation_pixel_sum[pixel_counter] += padded_reference.real_values[pixel_counter];
    //             //I am not sure about the slice to read value
    //             cc_output.WriteSlice(padded_reference.real_values[pixel_counter], pixel_counter);
    //             }
    //             padded_reference.SquareRealValues( );
    //             //                    correlation_pixel_sum_of_squares.AddImage(&padded_reference);
    //             for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated; pixel_counter++ ) {
    //                 correlation_pixel_sum_of_squares[pixel_counter] += padded_reference.real_values[pixel_counter];
    //                 }

    //                 //max_intensity_projection.QuickAndDirtyWriteSlice("/tmp/mip.mrc", 1);

    //                 current_projection.is_in_real_space = false;
    //                 padded_reference.is_in_real_space   = true;

    // }

    //what about code 890 - 1040?
    // write out one single MIP
    //is the mip just a single float value? MIP is an image max_intensity_projection object
    return true;
}
}