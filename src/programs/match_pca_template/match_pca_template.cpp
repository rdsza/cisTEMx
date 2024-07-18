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
Image test;

class
        MatchTemplateApp : public MyApp {
  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );
    void MasterHandleProgramDefinedResult(float* result_array, long array_size, int result_number, int number_of_expected_results);
    void ProgramSpecificInit( );

    // for master collation

    ArrayOfAggregatedTemplateResults aggregated_results;

    float GetMaxJobWaitTimeInSeconds( ) { return 360.0f; }

  private:
};

IMPLEMENT_APP(MatchTemplateApp)

void MatchTemplateApp::DoInteractiveUserInput( ) {
    UserInput* my_input = new UserInput("Perform Cross Correlations on image stack", 1.00);

    // get input
    wxString input_search_images      = my_input->GetFilenameFromUser("Input images to be searched", "", "input.mrc", false);
    wxString search_templates         = my_input->GetFilenameFromUser("Input template stack", "", "input.mrc", false);
    float defocus1                    = my_input->GetFloatFromUser("Defocus1 (angstroms)", "Defocus1 for the input image", "10000", 0.0);
    float defocus2                    = my_input->GetFloatFromUser("Defocus2 (angstroms)", "Defocus2 for the input image", "10000", 0.0);
    float defocus_angle               = my_input->GetFloatFromUser("Defocus Angle (degrees)", "Defocus Angle for the input image", "0.0");
    int first_search_position = 1; //is it one or zero
    int last_search_position = 91; //maybe 92
    //padding?
    //particle_radius_angstroms?
    delete my_input;

    my_current_job.ManualSetArguments("ttfffii", input_search_images.ToUTF8( ).data( ), search_templates.ToUTF8( ).data( ),
    defocus1, defocus2, defocus_angle, first_search_position, last_search_position);
  

}

bool MatchTemplateApp::DoCalculation( ) {
    //Bring inputs over from input function
    wxString input_search_images_filename  = my_current_job.arguments[0].ReturnStringArgument( );
    wxString input_reconstruction_filename = my_current_job.arguments[1].ReturnStringArgument( );
    float    defocus1                      = my_current_job.arguments[2].ReturnFloatArgument( );
    float    defocus2                      = my_current_job.arguments[3].ReturnFloatArgument( );
    float    defocus_angle                 = my_current_job.arguments[4].ReturnFloatArgument( );
    int      first_search_position         = my_current_job.arguments[5].ReturnIntegerArgument( );
    int      last_search_position          = my_current_job.arguments[6].ReturnIntegerArgument( );
    
    // do Template Match
    
    for ( int current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++ ) {
        // make the projection filter, which will be CTF * whitening filter
        input_ctf.SetDefocus((defocus1 + float(defocus_i) * defocus_step) / pixel_size, (defocus2 + float(defocus_i) * defocus_step) / pixel_size, deg_2_rad(defocus_angle));
        //            input_ctf.SetDefocus((defocus1 + 200) / pixel_size, (defocus2 + 200) / pixel_size, deg_2_rad(defocus_angle));
        projection_filter.CalculateCTFImage(input_ctf);
        projection_filter.ApplyCurveFilter(&whitening_filter);
        if ( padding != 1.0f ) {
                        template_reconstruction.ExtractSlice(padded_projection, angles, 1.0f, false); // TODO: change to ReadSlice
                        padded_projection.SwapRealSpaceQuadrants( );
                        padded_projection.BackwardFFT( );
                        padded_projection.ClipInto(&current_projection);
                        current_projection.ForwardFFT( );
                    }
                    else {
                        template_reconstruction.ExtractSlice(current_projection, angles, 1.0f, false); // TODO: change to ReadSlice
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
    // write out one single MIP
    
}