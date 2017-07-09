//
// The implementation is based on the following article:
// http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648
//

// x [-PI,PI]
float fast_sine(float x) {
    constexpr float PI = 3.14159265358f;
    constexpr float B = 4.0f / PI;
    constexpr float C = 4.0f / (PI * PI);
    constexpr float P = 0.225f;

    // always wrap input angle to -PI..PI
    if (x < -PI)
        x += 2*PI;
    else if (x > PI)
        x -= 2*PI;

    float y;
    if (x < 0) {
        y = B * x + C * x * x;

        if (y < 0)
            y = P * (y *-y - y) + y;
        else
            y = P * (y * y - y) + y;
    } else {
        y = B * x - C * x * x;

        if (y < 0)
            y = P * (y *-y - y) + y;
        else
            y = P * (y * y - y) + y;
    }
    return y;
}

// x [0, 2*PI]
float fast_cosine(float x) {
    constexpr float PI = 3.14159265358f;
    constexpr float B = 4.0f / PI;
    constexpr float C = 4.0f / (PI * PI);
    constexpr float P = 0.225f;

    //compute cosine: sin(x + PI/2) = cos(x)
    x += PI/2;
    if (x > PI)
        x -= 2*PI;

    float y;
    if (x < 0) {
        y = B * x + C * x * x;

        if (y < 0)
            y = P * (y *-y - y) + y;
        else
            y = P * (y * y - y) + y;
    } else {
        y = B * x - C * x * x;

        if (y < 0)
            y = P * (y *-y - y) + y;
        else
            y = P * (y * y - y) + y;
    }
    return y;
}

#include <cmath>
#include <cstdio>
#include <vector>

float calculate_std_deviation(const std::vector<float>& values, float* mean_value = nullptr) {
    float mean = 0.0;
    for (float v : values)
        mean += v;
    mean /= values.size();

    float variance = 0.0;
    for (float v : values)
        variance += (v - mean) * (v - mean);
    variance /= values.size();

    if (mean_value != nullptr)
        *mean_value = mean;

    return std::sqrt(variance);
}

int main() {
    const float PI = 3.14159265358f;
    const int N = 10'000;

    using Func = float(*)(float);

    auto calculate_errors = [&PI, &N](Func f_precise, Func f_fast, auto& abs_errors, auto& max_abs_error, auto& rel_errors, auto& max_rel_error) {
        max_abs_error = 0.0;
        max_rel_error = 0.0;
        for (int i = 0; i < N; i++) {
            float a = 2 * PI / (N - 1) * i;

            float precise = f_precise(a);
            float fast = f_fast(a);

            float abs_error = std::abs(precise - fast);
            abs_errors.push_back(abs_error);

            if (abs_error > max_abs_error)
                max_abs_error = abs_error;

            if (std::abs(precise) > 1e-5f) {
                float rel_error = abs_error / std::abs(precise);
                rel_errors.push_back(rel_error);

                if (rel_error > max_rel_error)
                    max_rel_error = rel_error;
            }
        }
    };

    //
    // Test sine accuracy.
    //
    {
        std::vector<float> abs_errors;
        float max_abs_error;
        std::vector<float> rel_errors;
        float max_rel_error;

        calculate_errors(std::sin, fast_sine, abs_errors, max_abs_error, rel_errors, max_rel_error);

        float abs_error_mean;
        float abs_error_std_deviation = calculate_std_deviation(abs_errors, &abs_error_mean);
        float rel_error_mean;
        float rel_error_std_deviation = calculate_std_deviation(rel_errors, &rel_error_mean);

        printf("sine max abs error = %f\n", max_abs_error);
        printf("sine abs error mean = %f\n", abs_error_mean);
        printf("sine abs error std deviation = %f\n\n", abs_error_std_deviation);

        printf("sine max relative error = %.3f%%\n", max_rel_error * 100.0f);
        printf("sine relative error mean = %.3f%%\n", rel_error_mean * 100.0f);
        printf("sine relative error std deviation = %.3f%%\n\n", rel_error_std_deviation * 100.0f);
    }

    //
    // Test cosine accuracy.
    //
    {
        std::vector<float> abs_errors;
        float max_abs_error;
        std::vector<float> rel_errors;
        float max_rel_error;

        calculate_errors(std::cos, fast_cosine, abs_errors, max_abs_error, rel_errors, max_rel_error);

        float abs_error_mean;
        float abs_error_std_deviation = calculate_std_deviation(abs_errors, &abs_error_mean);
        float rel_error_mean;
        float rel_error_std_deviation = calculate_std_deviation(rel_errors, &rel_error_mean);

        printf("cosine max abs error = %f\n", max_abs_error);
        printf("cosine abs error mean = %f\n", abs_error_mean);
        printf("cosine abs error std deviation = %f\n\n", abs_error_std_deviation);

        printf("cosine max relative error = %.3f%%\n", max_rel_error * 100.0f);
        printf("cosine relative error mean = %.3f%%\n", rel_error_mean * 100.0f);
        printf("cosine relative error std deviation = %.3f%%\n", rel_error_std_deviation * 100.0f);
    }
}
