//
// The implementation is based on the following article:
// http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648
//

// x [-PI,PI]
float fast_sine(float x) {
    constexpr float PI = 3.14159265358f;
    constexpr float B = 4.0f / PI;
    constexpr float C = -4.0f / (PI * PI);
    constexpr float P = 0.225f;

    float y = B * x + C * x * (x < 0 ? -x : x);
    return P * (y * (y < 0 ? -y : y) - y) + y;
}

// x [0, 2*PI]
float fast_cosine(float x) {
    constexpr float PI = 3.14159265358f;
    constexpr float B = 4.0f / PI;
    constexpr float C = -4.0f / (PI * PI);
    constexpr float P = 0.225f;

    x += PI/2;
    if (x > PI)
        x -= 2*PI;

    return fast_sine(x);
}

#include <cmath>
#include <chrono>
#include <cstdio>
#include <vector>

int main() {
    auto calculate_std_deviation = [](const std::vector<float>& values, float* mean_value = nullptr) {
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
    };

    //
    // Test sine accuracy.
    //
    {
        const float PI = 3.14159265358f;
        const int N = 10'000;

        std::vector<float> abs_errors;
        float max_abs_error = 0.0f;
        std::vector<float> rel_errors;
        float max_rel_error = 0.0f;

        for (int i = 0; i < N; i++) {
            float a = -PI + 2 * PI / (N - 1) * i;

            float precise = std::sin(a);
            float fast = fast_sine(a);

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
        const float PI = 3.14159265358f;
        const int N = 10'000;

        std::vector<float> abs_errors;
        float max_abs_error = 0.0f;
        std::vector<float> rel_errors;
        float max_rel_error = 0.0f;

        for (int i = 0; i < N; i++) {
            float a = 2 * PI / (N - 1) * i;

            float precise = std::cos(a);
            float fast = fast_cosine(a);

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

        float abs_error_mean;
        float abs_error_std_deviation = calculate_std_deviation(abs_errors, &abs_error_mean);
        float rel_error_mean;
        float rel_error_std_deviation = calculate_std_deviation(rel_errors, &rel_error_mean);

        printf("cosine max abs error = %f\n", max_abs_error);
        printf("cosine abs error mean = %f\n", abs_error_mean);
        printf("cosine abs error std deviation = %f\n\n", abs_error_std_deviation);

        printf("cosine max relative error = %.3f%%\n", max_rel_error * 100.0f);
        printf("cosine relative error mean = %.3f%%\n", rel_error_mean * 100.0f);
        printf("cosine relative error std deviation = %.3f%%\n\n", rel_error_std_deviation * 100.0f);
    }

    //
    // Test performance.
    //
    {
        struct Timer {
            using Clock = std::chrono::high_resolution_clock;
            Clock::time_point start = Clock::now();
            int elapsed_microseconds() const {
                const auto duration = Clock::now() - start;
                return (int)std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
            }
        };

        const float PI = 3.14159265358f;
        const int N = 1'000'000;
        const float da = 2*PI / (N - 1);

        // sine
        {
            Timer t;
            float s = 0, a = -PI;
            for (int i = 0; i < N; i++, a += da) {
                s += std::sin(a);
            }
            auto elapsed = t.elapsed_microseconds();
            printf("Regular sine time: %d\n", elapsed);

            if (s == PI) printf("Miracle"); // this prevents the compiler from optimizing out the test loop
        }

        {
            Timer t;
            float s = 0, a = -PI;
            for (int i = 0; i < N; i++, a += da) {
                s += fast_sine(a);
            }
            auto elapsed = t.elapsed_microseconds();
            printf("Fast sine time: %d\n\n", elapsed);

            if (s == PI) printf("Miracle"); 
        }

        // cosine
        {
            Timer t;
            float s = 0, a = 0;
            for (int i = 0; i < N; i++, a += da) {
                s += std::cos(a);
            }
            auto elapsed = t.elapsed_microseconds();
            printf("Regular cosine time: %d\n", elapsed);

            if (s == PI) printf("Miracle");
        }

        {
            Timer t;
            float s = 0, a = 0;
            for (int i = 0; i < N; i++, a += da) {
                s += fast_cosine(a);
            }
            auto elapsed = t.elapsed_microseconds();
            printf("Fast cosine time: %d\n", elapsed);

            if (s == PI) printf("Miracle");
        }
    }
}
