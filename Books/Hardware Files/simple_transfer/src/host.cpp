////
////// Basic Host Code to run one kernel
////#include <iostream>
////#include <iomanip>
////#include <cstring>
////#include <string>
////
////// XRT includes
////#include "experimental/xrt_device.h"
////#include "experimental/xrt_kernel.h"
////#include "experimental/xrt_bo.h"
////
////#define MAX_SEQ_LEN 128
////#define ALIGN_LEN (MAX_SEQ_LEN * 2)
////
////// Helper: Trim trailing nulls or spaces
////std::string trim_nulls(const char* str, size_t len) {
////    size_t end = len;
////    while (end > 0 && (str[end - 1] == '\0' || str[end - 1] == ' '))
////        --end;
////    return std::string(str, end);
////}
////
////int main(int argc, char** argv) {
////    // Load basic arguments
////    if (argc < 2) {
////        std::cerr << "Usage: " << argv[0] << " <xclbin file>\n";
////        return EXIT_FAILURE;
////    }
////
////    std::string xclbin_file = argv[1];
////
////    // User Input
////    std::string seq1, seq2;
////    std::cout << "🔤 Enter DNA Sequence 1 (max " << MAX_SEQ_LEN << " chars): ";
////    std::getline(std::cin, seq1);
////    std::cout << "🔤 Enter DNA Sequence 2 (max " << MAX_SEQ_LEN << " chars): ";
////    std::getline(std::cin, seq2);
////
////    if (seq1.length() > MAX_SEQ_LEN || seq2.length() > MAX_SEQ_LEN) {
////        std::cerr << "❌ Error: Input sequences exceed max length of " << MAX_SEQ_LEN << " characters.\n";
////        return EXIT_FAILURE;
////    }
////
////    int len1 = seq1.length();
////    int len2 = seq2.length();
////
////    int gap_penalty = -2;
////    int gap_extend = -1;
////    int match_score = 3;
////    int mismatch_penalty = -3;
////
////    std::cout << "📡 Connecting to device...\n";
////    auto device = xrt::device(0);
////    auto uuid = device.load_xclbin(xclbin_file);
////    auto kernel = xrt::kernel(device, uuid, "bwfa_kernel");
////
////    // Allocate buffers
////    auto bo_seq1 = xrt::bo(device, MAX_SEQ_LEN, kernel.group_id(0));
////    auto bo_seq2 = xrt::bo(device, MAX_SEQ_LEN, kernel.group_id(1));
////    auto bo_aligned1 = xrt::bo(device, ALIGN_LEN, kernel.group_id(8));
////    auto bo_aligned2 = xrt::bo(device, ALIGN_LEN, kernel.group_id(9));
////
////    char* ptr_seq1 = bo_seq1.map<char*>();
////    char* ptr_seq2 = bo_seq2.map<char*>();
////    char* ptr_aligned1 = bo_aligned1.map<char*>();
////    char* ptr_aligned2 = bo_aligned2.map<char*>();
////
////    std::memset(ptr_seq1, 0, MAX_SEQ_LEN);
////    std::memset(ptr_seq2, 0, MAX_SEQ_LEN);
////    std::memset(ptr_aligned1, 0, ALIGN_LEN);
////    std::memset(ptr_aligned2, 0, ALIGN_LEN);
////
////    std::memcpy(ptr_seq1, seq1.c_str(), len1);
////    std::memcpy(ptr_seq2, seq2.c_str(), len2);
////
////    bo_seq1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
////    bo_seq2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
////
////    std::cout << "🚀 Launching kernel...\n";
////    auto run = kernel(
////        bo_seq1, bo_seq2,
////        len1, len2,
////        gap_penalty, gap_extend, match_score, mismatch_penalty,
////        bo_aligned1, bo_aligned2
////    );
////    run.wait();
////
////    bo_aligned1.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
////    bo_aligned2.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
////
////    std::string aligned1 = trim_nulls(ptr_aligned1, ALIGN_LEN);
////    std::string aligned2 = trim_nulls(ptr_aligned2, ALIGN_LEN);
////
////    std::cout << "\n==================[ ALIGNMENT RESULT ]==================\n";
////    std::cout << "🔬 Input Seq1  : " << seq1 << "\n";
////    std::cout << "🔬 Input Seq2  : " << seq2 << "\n";
////    std::cout << "--------------------------------------------------------\n";
////    std::cout << "✅ Aligned 1   : " << aligned1 << "\n";
////    std::cout << "✅ Aligned 2   : " << aligned2 << "\n";
////    std::cout << "========================================================\n";
////
////    return 0;
////}
//
//// Affine- Gap Host Code to run 4 kernels parallely
//#include <iostream>
//#include <iomanip>
//#include <cstring>
//#include <string>
//#include <vector>
//#include <memory>
//#include <algorithm>
//#include <time.h>
//
//// XRT includes
//#include "experimental/xrt_device.h"
//#include "experimental/xrt_kernel.h"
//#include "experimental/xrt_bo.h"
//
//#define MAX_SEQ_LEN 128
//#define ALIGN_LEN (MAX_SEQ_LEN * 2)
//#define NUM_KERNELS 4  // Number of CUs in xclbin
//
//// Helper: Trim trailing nulls or spaces
//std::string trim_nulls(const char* str, size_t len) {
//    size_t end = len;
//    while (end > 0 && (str[end - 1] == '\0' || str[end - 1] == ' '))
//        --end;
//    return std::string(str, end);
//}
//
////gettimestamp: Function to get the execution time of the kernel
//double gettimestamp(){
//	struct timeval tv;
//	gettimeofday(&tv,NULL);
//	return tv.tv_usec + tv.tv_sec*1e6;
//}
//
//double hardware_start;
//double hardware_end;
//double hardware_time;
//
//int main(int argc, char** argv) {
//    if (argc < 2) {
//        std::cerr << "Usage: " << argv[0] << " <xclbin file>\n";
//        return EXIT_FAILURE;
//    }
//
//    std::string xclbin_file = argv[1];
//
//    // User Input
//    std::string seq1, seq2;
//    std::cout << "🔤 Enter DNA Sequence 1 (max " << MAX_SEQ_LEN << " chars): ";
//    std::getline(std::cin, seq1);
//    std::cout << "🔤 Enter DNA Sequence 2 (max " << MAX_SEQ_LEN << " chars): ";
//    std::getline(std::cin, seq2);
//
//    if (seq1.length() > NUM_KERNELS * MAX_SEQ_LEN || seq2.length() > NUM_KERNELS * MAX_SEQ_LEN) {
//        std::cerr << "❌ Error: Input sequences exceed max length of " << NUM_KERNELS * MAX_SEQ_LEN << " characters.\n";
//        return EXIT_FAILURE;
//    }
//
//    int len1 = seq1.length();
//    int len2 = seq2.length();
//
//    int gap_penalty = -2;
//    int gap_extend = -1;
//    int match_score = 1;
//    int mismatch_penalty = -2;
//
//    std::cout << "📡 Connecting to device...\n";
//    auto device = xrt::device(0);
//
//    std::cout << "Starting Timer\n";
//    hardware_start = gettimestamp();
//
//    auto uuid = device.load_xclbin(xclbin_file);
//
//    //int part_len1 = (len1 + NUM_KERNELS - 1) / NUM_KERNELS;
//    //int part_len2 = (len2 + NUM_KERNELS - 1) / NUM_KERNELS;
//
//    const int part_len = MAX_SEQ_LEN;  //Fixed 128 chars per kernel
//
//
//    std::vector<xrt::bo> bos_seq1, bos_seq2, bos_aligned1, bos_aligned2;
//    std::vector<std::unique_ptr<char[]>> maps_seq1, maps_seq2;
//    std::vector<xrt::run> runs;
//
//    for (int i = 0; i < NUM_KERNELS; ++i) {
//
//    	// Calculate part lengths and start positions
//    	int start1 = i * part_len;
//    	int len_chunk1 = std::min(part_len, (len1 > start1) ? (len1 - start1) : 0 );
//    	int start2 = i * part_len;
//    	int len_chunk2 = std::min(part_len, (len2 > start2) ? (len2 - start2) : 0 );
//
//    	// Checking inputs for the kernels
//    	if(len_chunk1 == 0 || len_chunk2 == 0){
//    		std::cout << "Skipping kernel instance " << i + 1 << " due to zero length input.\n";
//    		continue;
//    	}
//
//        // Use CU-specific kernel name
//        std::string kernel_name = "bwfa_kernel:{bwfa_kernel_" + std::to_string(i + 1) + "}";
//        auto kernel = xrt::kernel(device, uuid, kernel_name);
//
//        // Allocate buffers
//        bos_seq1.emplace_back(device, MAX_SEQ_LEN, kernel.group_id(0));
//        bos_seq2.emplace_back(device, MAX_SEQ_LEN, kernel.group_id(1));
//        bos_aligned1.emplace_back(device, ALIGN_LEN, kernel.group_id(8));
//        bos_aligned2.emplace_back(device, ALIGN_LEN, kernel.group_id(9));
//
//        maps_seq1.emplace_back(new char[MAX_SEQ_LEN]());
//        maps_seq2.emplace_back(new char[MAX_SEQ_LEN]());
//
//        std::memcpy(maps_seq1[i].get(), seq1.c_str() + start1, len_chunk1);
//        std::memcpy(maps_seq2[i].get(), seq2.c_str() + start2, len_chunk2);
//
//        std::memcpy(bos_seq1[i].map<char*>(), maps_seq1[i].get(), MAX_SEQ_LEN);
//        std::memcpy(bos_seq2[i].map<char*>(), maps_seq2[i].get(), MAX_SEQ_LEN);
//
//        bos_seq1[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
//        bos_seq2[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
//
//        // Launch the kernel
//        std::cout << "🚀 Launching kernel instance " << i + 1 << "...\n";
//        runs.emplace_back(kernel(
//            bos_seq1[i], bos_seq2[i],
//            len_chunk1, len_chunk2,
//            gap_penalty, gap_extend, match_score, mismatch_penalty,
//            bos_aligned1[i], bos_aligned2[i]
//        ));
//    }
//
//    // Wait for all runs to finish
//    for (auto& run : runs)
//        run.wait();
//
//    std::cout << "Ending Timer\n";
//    hardware_end = gettimestamp();
//
//    hardware_time = (hardware_end - hardware_start)/1000;
//
//    // Retrieve aligned results
//    std::string final_aligned1, final_aligned2;
//    for (int i = 0; i < runs.size(); ++i) {
//        bos_aligned1[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
//        bos_aligned2[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
//
//        char* p1 = bos_aligned1[i].map<char*>();
//        char* p2 = bos_aligned2[i].map<char*>();
//
//        final_aligned1 += trim_nulls(p1, ALIGN_LEN);
//        final_aligned2 += trim_nulls(p2, ALIGN_LEN);
//    }
//
//    std::cout << "\n==================[ ALIGNMENT RESULT ]==================\n";
//    std::cout << "🔬 Input Seq1  : " << seq1 << "\n";
//    std::cout << "🔬 Input Seq2  : " << seq2 << "\n";
//    std::cout << "--------------------------------------------------------\n";
//    std::cout << "✅ Aligned 1   : " << final_aligned1 << "\n";
//    std::cout << "✅ Aligned 2   : " << final_aligned2 << "\n";
//    std::cout << "========================================================\n";
//    std::cout << "Execution Time: " << hardware_time << "msec\n";
//    std::cout << "========================================================\n";
//
//    return 0;
//}
//
//
//
//

//AFFINE Gap With K-Mer Logic

//#include <iostream>
//#include <iomanip>
//#include <cstring>
//#include <string>
//#include <vector>
//#include <memory>
//#include <algorithm>
//#include <time.h>
//#include <fstream>
//
//// XRT includes
//#include "experimental/xrt_device.h"
//#include "experimental/xrt_kernel.h"
//#include "experimental/xrt_bo.h"
//
//#define MAX_SEQ_LEN 64
//#define ALIGN_LEN (MAX_SEQ_LEN * 2)
//#define TOTAL_KERNELS 4
//#define KMER_SIZE 16
//
//std::string trim_nulls(const char* str, size_t len) {
//    size_t end = len;
//    while (end > 0 && (str[end - 1] == '\0' || str[end - 1] == ' '))
//        --end;
//    return std::string(str, end);
//}
//
//std::vector<std::string> generate_kmers(const std::string& sequence, int k) {
//    std::vector<std::string> kmers;
//    for (size_t i = 0; i <= sequence.length() - k; ++i) {
//        kmers.push_back(sequence.substr(i, k));
//    }
//    return kmers;
//}
//
//std::pair<std::string, std::string> merge_with_best_kmer_overlap(
//    const std::vector<std::string>& chunks1,
//    const std::vector<std::string>& chunks2,
//    size_t min_kmer_overlap = 8
//) {
//    if (chunks1.empty()) return {"", ""};
//    std::string result1 = chunks1[0];
//    std::string result2 = chunks2[0];
//
//    for (size_t i = 1; i < chunks1.size(); ++i) {
//        const std::string& next1 = chunks1[i];
//        const std::string& next2 = chunks2[i];
//
//        size_t max_overlap = std::min(result1.size(), next1.size());
//        size_t best_overlap = 0;
//
//        // Find the longest overlap from end of result1/2 and start of next1/2
//        for (size_t k = max_overlap; k >= min_kmer_overlap; --k) {
//            if (
//                result1.substr(result1.size() - k) == next1.substr(0, k) &&
//                result2.substr(result2.size() - k) == next2.substr(0, k)
//            ) {
//                best_overlap = k;
//                break;
//            }
//        }
//
//        // Append non-overlapping parts only
//        result1 += next1.substr(best_overlap);
//        result2 += next2.substr(best_overlap);
//    }
//
//    return {result1, result2};
//}
//
//
//bool read_sequence_from_file(const std::string& filename, std::string& sequence) {
//    std::ifstream file(filename);
//    if (!file.is_open()) {
//        std::cerr << "❌ Error: Could not open the file " << filename << "\n";
//        return false;
//    }
//    std::getline(file, sequence);
//    file.close();
//    if (sequence.empty()) {
//        std::cerr << "❌ Error: The sequence in file " << filename << " is empty.\n";
//        return false;
//    }
//    return true;
//}
//
//double gettimestamp() {
//    struct timeval tv;
//    gettimeofday(&tv, NULL);
//    return tv.tv_usec + tv.tv_sec * 1e6;
//}
////double hardware_start;
////double hardware_end;
////double hardware_time;
////double hardware_kerneltime;
//double hardware_kerneltotal=0;
//double hardware_kernelmax = 0;

//int main(int argc, char** argv) {
//    if (argc < 2) {
//        std::cerr << "Usage: " << argv[0] << " <xclbin file>\n";
//        return EXIT_FAILURE;
//    }
//
//    std::string xclbin_file = argv[1];
//    std::string seq_file1, seq_file2;
//    std::cout << "🔤 Enter the address of the file for DNA Sequence 1: ";
//    std::getline(std::cin, seq_file1);
//    std::cout << "🔤 Enter the address of the file for DNA Sequence 2: ";
//    std::getline(std::cin, seq_file2);
//
//    std::string seq1, seq2;
//    if (!read_sequence_from_file(seq_file1, seq1)) return EXIT_FAILURE;
//    if (!read_sequence_from_file(seq_file2, seq2)) return EXIT_FAILURE;
//
//    int len1 = seq1.length();
//    int len2 = seq2.length();
//    int gap_penalty = -2, gap_extend = -1, match_score = 1, mismatch_penalty = -2;
//
//    if (len1 == 0 || len2 == 0) {
//        std::cerr << "❌ Error: One or both sequences are empty.\n";
//        return EXIT_FAILURE;
//    }
//
//    auto kmers1 = generate_kmers(seq1, KMER_SIZE);
//    auto kmers2 = generate_kmers(seq2, KMER_SIZE);
//
//    int num_kmers = std::min(kmers1.size(), kmers2.size());
//    int total_kernels = 4;
//    int num_iterations = (num_kmers + total_kernels - 1) / total_kernels;
//
//    std::string final_aligned1, final_aligned2;
//
//    std::cout << "\n====================[ DNA SEQUENCE ALIGNMENT ]====================\n";
//    std::cout << "🔑 Sequence 1 Length: " << len1 << " | Sequence 2 Length: " << len2 << "\n";
//    std::cout << "⚙️  Gap Penalty: " << gap_penalty << "Gap Extend: " << gap_extend <<", Match Score: " << match_score << ", Mismatch Penalty: " << mismatch_penalty << "\n";
//    std::cout << "===================================================================\n";
//
//    std::cout << "📡 Connecting to FPGA device...\n";
//    auto device = xrt::device(0);
//    std::cout << "⚡ Device Name: " << device.get_info<xrt::info::device::name>() << "\n";
//
//    auto uuid = device.load_xclbin(xclbin_file);
//
//    // Allocate buffers statically
//    std::vector<xrt::bo> bos_seq1(TOTAL_KERNELS), bos_seq2(TOTAL_KERNELS), bos_aligned1(TOTAL_KERNELS), bos_aligned2(TOTAL_KERNELS);
//    std::vector<std::unique_ptr<char[]>> maps_seq1(TOTAL_KERNELS), maps_seq2(TOTAL_KERNELS);
//    std::vector<xrt::kernel> kernels(TOTAL_KERNELS);
//
//    for (int i = 0; i < TOTAL_KERNELS; ++i) {
//        std::string kernel_name = "bwfa_kernel:{bwfa_kernel_" + std::to_string(i + 1) + "}";
//        kernels[i] = xrt::kernel(device, uuid, kernel_name);
//
//        bos_seq1[i] = xrt::bo(device, MAX_SEQ_LEN, kernels[i].group_id(0));
//        bos_seq2[i] = xrt::bo(device, MAX_SEQ_LEN, kernels[i].group_id(1));
//        bos_aligned1[i] = xrt::bo(device, ALIGN_LEN, kernels[i].group_id(8));
//        bos_aligned2[i] = xrt::bo(device, ALIGN_LEN, kernels[i].group_id(9));
//
//        maps_seq1[i] = std::make_unique<char[]>(MAX_SEQ_LEN);
//        maps_seq2[i] = std::make_unique<char[]>(MAX_SEQ_LEN);
//    }
//
//    std::cout << "Starting Timer: " << "\n";
//    double hardware_start = gettimestamp();
//    double hardware_kerneltime = 0;
//
//
//    for (int iteration = 0; iteration < num_iterations; ++iteration) {
//    std::cout << "Starting Iteration: " << iteration + 1 << "\n";
//    runs.clear();
//    bos_seq1.clear(); bos_seq2.clear();
//    bos_aligned1.clear(); bos_aligned2.clear();
//    maps_seq1.clear(); maps_seq2.clear();
//
//    for (int i = 0; i < total_kernels; ++i) {
//        int kmer_index = iteration * total_kernels + i;
//        if (kmer_index >= num_kmers) break;
//
//        std::string kmer1 = kmers1[kmer_index];
//        std::string kmer2 = kmers2[kmer_index];
//
//        int len_kmer1 = kmer1.length();
//        int len_kmer2 = kmer2.length();
//
//        //std::string kernel_name = "bwfa_kernel:{bwfa_kernel_" + std::to_string(i + 1) + "}";
//        //auto kernel = xrt::kernel(device, uuid, kernel_name);
//
//        bos_seq1.emplace_back(device, MAX_SEQ_LEN, kernel.group_id(0));
//        bos_seq2.emplace_back(device, MAX_SEQ_LEN, kernel.group_id(1));
//        bos_aligned1.emplace_back(device, ALIGN_LEN, kernel.group_id(8));
//        bos_aligned2.emplace_back(device, ALIGN_LEN, kernel.group_id(9));
//
//        maps_seq1.emplace_back(new char[MAX_SEQ_LEN]());
//        maps_seq2.emplace_back(new char[MAX_SEQ_LEN]());
//
//        std::memcpy(maps_seq1[i].get(), kmer1.c_str(), len_kmer1);
//        std::memcpy(maps_seq2[i].get(), kmer2.c_str(), len_kmer2);
//
//        std::memcpy(bos_seq1[i].map<char*>(), maps_seq1[i].get(), MAX_SEQ_LEN);
//        std::memcpy(bos_seq2[i].map<char*>(), maps_seq2[i].get(), MAX_SEQ_LEN);
//
//        bos_seq1[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
//        bos_seq2[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
//
//        std::cout << "🚀 Launching kernel instance " << i + 1 << "...\n";
//
//        double t_start = gettimestamp();
//        runs.emplace_back(kernel(
//            bos_seq1[i], bos_seq2[i],
//            len_kmer1, len_kmer2,
//            gap_penalty, gap_extend, match_score, mismatch_penalty,
//            bos_aligned1[i], bos_aligned2[i]
//        ));
//        double t_end = gettimestamp();
//        hardware_kerneltime = (t_end - t_start) / 1000.0;
//        hardware_kernelmax = (hardware_kernelmax > hardware_kerneltime) ? hardware_kernelmax : hardware_kerneltime;
//        std::cout << "Kernel Execution Time: "<< hardware_kerneltime << "msec\n";
//        runs.emplace_back(std::move(run));
//    }
//
//    for (auto& run : runs) run.wait();
//
//    std::vector<std::string> aligned_chunks1;
//    std::vector<std::string> aligned_chunks2;
//
//    for (int i = 0; i < runs.size(); ++i) {
//        bos_aligned1[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
//        bos_aligned2[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
//
//        char* p1 = bos_aligned1[i].map<char*>();
//        char* p2 = bos_aligned2[i].map<char*>();
//
//        aligned_chunks1.push_back(trim_nulls(p1, ALIGN_LEN));
//        aligned_chunks2.push_back(trim_nulls(p2, ALIGN_LEN));
//    }
//
//    // Merge using overlap-aware logic
//    std::tie(final_aligned1, final_aligned2) = merge_with_best_kmer_overlap(aligned_chunks1, aligned_chunks2);
//
//    hardware_kerneltotal=hardware_kerneltotal+hardware_kernelmax;
//}
//
//    double hardware_end = gettimestamp();
//    double total_time = (hardware_end - hardware_start) / 1000;
//    std::cout << "Ending timer" << "\n";
//
//    std::cout << "\n==================[ ALIGNMENT RESULT ]==================\n";
//    std::cout << "🔬 Input Seq1  : " << seq1 << "\n";
//    std::cout << "🔬 Input Seq2  : " << seq2 << "\n";
//    std::cout << "--------------------------------------------------------\n";
//    std::cout << "✅ Aligned 1   : " << final_aligned1 << "\n";
//    std::cout << "✅ Aligned 2   : " << final_aligned2 << "\n";
//    std::cout << "========================================================\n";
//    std::cout << "🕒 Total Execution Time: " << total_time << " msec\n";
//    std::cout << "⚙️  Total Kernel-only Time   : " << hardware_kerneltotal << " msec\n";
//    std::cout << "========================================================\n";
//    std::cout << "🎉 Alignment completed successfully!\n";
//
//    return 0;
//}



#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <time.h>
#include <fstream>

// XRT includes
#include "experimental/xrt_device.h"
#include "experimental/xrt_kernel.h"
#include "experimental/xrt_bo.h"

#define MAX_SEQ_LEN 64
#define ALIGN_LEN (MAX_SEQ_LEN * 2)
#define TOTAL_KERNELS 4

std::string trim_nulls(const char* str, size_t len) {
    size_t end = len;
    while (end > 0 && (str[end - 1] == '\0' || str[end - 1] == ' '))
        --end;
    return std::string(str, end);
}

bool read_sequence_from_file(const std::string& filename, std::string& sequence) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "❌ Error: Could not open the file " << filename << "\n";
        return false;
    }
    std::getline(file, sequence);
    file.close();
    if (sequence.empty()) {
        std::cerr << "❌ Error: The sequence in file " << filename << " is empty.\n";
        return false;
    }
    return true;
}

double gettimestamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_usec + tv.tv_sec * 1e6;
}
//double hardware_start;
//double hardware_end;
//double hardware_time;
//double hardware_kerneltime;
double hardware_kerneltotal=0;
double hardware_kernelmax = 0;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <xclbin file>\n";
        return EXIT_FAILURE;
    }

    std::string xclbin_file = argv[1];
    std::string seq_file1, seq_file2;
    std::cout << "🔤 Enter the address of the file for DNA Sequence 1: ";
    std::getline(std::cin, seq_file1);
    std::cout << "🔤 Enter the address of the file for DNA Sequence 2: ";
    std::getline(std::cin, seq_file2);

    std::string seq1, seq2;
    if (!read_sequence_from_file(seq_file1, seq1)) return EXIT_FAILURE;
    if (!read_sequence_from_file(seq_file2, seq2)) return EXIT_FAILURE;

    int len1 = seq1.length();
    int len2 = seq2.length();
    int gap_penalty = -2, gap_extend = -1, match_score = 1, mismatch_penalty = -2;

    if (len1 == 0 || len2 == 0) {
        std::cerr << "❌ Error: One or both sequences are empty.\n";
        return EXIT_FAILURE;
    }

    int total_length = std::max(len1, len2);
    int num_iterations = (total_length + (TOTAL_KERNELS * MAX_SEQ_LEN - 1)) / (TOTAL_KERNELS * MAX_SEQ_LEN);
    std::string final_aligned1, final_aligned2;

    std::cout << "\n====================[ DNA SEQUENCE ALIGNMENT ]====================\n";
    std::cout << "🔑 Sequence 1 Length: " << len1 << " | Sequence 2 Length: " << len2 << "\n";
    std::cout << "⚙️  Gap Penalty: " << gap_penalty << ", Gap Extend: " << gap_extend <<", Match Score: " << match_score << ", Mismatch Penalty: " << mismatch_penalty << "\n";
    std::cout << "===================================================================\n";

    std::cout << "📡 Connecting to FPGA device...\n";
    auto device = xrt::device(0);
    std::cout << "⚡ Device Name: " << device.get_info<xrt::info::device::name>() << "\n";

    auto uuid = device.load_xclbin(xclbin_file);

    // Allocate buffers statically
    std::vector<xrt::bo> bos_seq1(TOTAL_KERNELS), bos_seq2(TOTAL_KERNELS), bos_aligned1(TOTAL_KERNELS), bos_aligned2(TOTAL_KERNELS);
    std::vector<std::unique_ptr<char[]>> maps_seq1(TOTAL_KERNELS), maps_seq2(TOTAL_KERNELS);
    std::vector<xrt::kernel> kernels(TOTAL_KERNELS);

    for (int i = 0; i < TOTAL_KERNELS; ++i) {
        std::string kernel_name = "bwfa_kernel:{bwfa_kernel_" + std::to_string(i + 1) + "}";
        kernels[i] = xrt::kernel(device, uuid, kernel_name);

        bos_seq1[i] = xrt::bo(device, MAX_SEQ_LEN, kernels[i].group_id(0));
        bos_seq2[i] = xrt::bo(device, MAX_SEQ_LEN, kernels[i].group_id(1));
        bos_aligned1[i] = xrt::bo(device, ALIGN_LEN, kernels[i].group_id(8));
        bos_aligned2[i] = xrt::bo(device, ALIGN_LEN, kernels[i].group_id(9));

        maps_seq1[i] = std::make_unique<char[]>(MAX_SEQ_LEN);
        maps_seq2[i] = std::make_unique<char[]>(MAX_SEQ_LEN);
    }

    std::cout << "Starting Timer: " << "\n";
    double hardware_start = gettimestamp();
    double hardware_kerneltime = 0;


    for (int iteration = 0; iteration < num_iterations; ++iteration) {
        std::cout << "▶ Iteration " << iteration + 1 << "\n";
        int start_offset = iteration * TOTAL_KERNELS * MAX_SEQ_LEN;
        std::vector<xrt::run> runs;
        runs.clear();
        hardware_kernelmax = 0;


        for (int i = 0; i < TOTAL_KERNELS; ++i) {
            int start1 = start_offset + i * MAX_SEQ_LEN;
            int start2 = start_offset + i * MAX_SEQ_LEN;

            int len_chunk1 = std::min(MAX_SEQ_LEN, (len1 > start1) ? (len1 - start1) : 0);
            int len_chunk2 = std::min(MAX_SEQ_LEN, (len2 > start2) ? (len2 - start2) : 0);

            if (len_chunk1 == 0 || len_chunk2 == 0) {
            	std::cout << "⏩ Skipping kernel instance " << i + 1 << " due to zero length input.\n";
                continue;
            }

            // Clear and copy input buffers
            std::memset(maps_seq1[i].get(), 0, MAX_SEQ_LEN);
            std::memset(maps_seq2[i].get(), 0, MAX_SEQ_LEN);
            std::memcpy(maps_seq1[i].get(), seq1.c_str() + start1, len_chunk1);
            std::memcpy(maps_seq2[i].get(), seq2.c_str() + start2, len_chunk2);

            std::memcpy(bos_seq1[i].map<char*>(), maps_seq1[i].get(), MAX_SEQ_LEN);
            std::memcpy(bos_seq2[i].map<char*>(), maps_seq2[i].get(), MAX_SEQ_LEN);

            // Clear output buffers to prevent stale data
            std::memset(bos_aligned1[i].map<char*>(), 0, ALIGN_LEN);
            std::memset(bos_aligned2[i].map<char*>(), 0, ALIGN_LEN);

            bos_seq1[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
            bos_seq2[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);

            bos_aligned1[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
            bos_aligned2[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
            // Launch kernel for current chunk of sequences
            std::cout << "🚀 Launching kernel instance " << i + 1 << "...\n";

            double t_start = gettimestamp();
            auto run = kernels[i](
                bos_seq1[i], bos_seq2[i],
                len_chunk1, len_chunk2,
                gap_penalty, gap_extend, match_score, mismatch_penalty,
                bos_aligned1[i], bos_aligned2[i]
            );
            double t_end = gettimestamp();
            hardware_kerneltime = (t_end - t_start) / 1000.0;
            hardware_kernelmax = (hardware_kernelmax > hardware_kerneltime) ? hardware_kernelmax : hardware_kerneltime;
            std::cout << "Kernel Execution Time: "<< hardware_kerneltime << "msec\n";
            runs.emplace_back(std::move(run));
        }

        std::cout << "Max Kernel Time For Iteration "<< iteration + 1 <<": "<< hardware_kernelmax << "msec\n";

        for (auto& run : runs) {
            run.wait();
        }

        for (int i = 0; i < TOTAL_KERNELS; ++i) {
            int start1 = start_offset + i * MAX_SEQ_LEN;
            int len_chunk1 = std::min(MAX_SEQ_LEN, (len1 > start1) ? (len1 - start1) : 0);
            if (len_chunk1 == 0) continue;

            bos_aligned1[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            bos_aligned2[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            final_aligned1 += trim_nulls(bos_aligned1[i].map<char*>(), ALIGN_LEN);
            final_aligned2 += trim_nulls(bos_aligned2[i].map<char*>(), ALIGN_LEN);
        }

        hardware_kerneltotal=hardware_kerneltotal+hardware_kernelmax;
    }

    double hardware_end = gettimestamp();
    double total_time = (hardware_end - hardware_start) / 1000;
    std::cout << "Ending timer" << "\n";

    std::cout << "\n==================[ ALIGNMENT RESULT ]==================\n";
    std::cout << "🔬 Input Seq1  : " << seq1 << "\n";
    std::cout << "🔬 Input Seq2  : " << seq2 << "\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << "✅ Aligned 1   : " << final_aligned1 << "\n";
    std::cout << "✅ Aligned 2   : " << final_aligned2 << "\n";
    std::cout << "========================================================\n";
    std::cout << "🕒 Total Execution Time: " << total_time << " msec\n";
//    std::cout << "⚙️  Total Kernel-only Time   : " << hardware_kerneltotal << " msec\n";
    std::cout << "========================================================\n";
    std::cout << "🎉 Alignment completed successfully!\n";

    return 0;
}
