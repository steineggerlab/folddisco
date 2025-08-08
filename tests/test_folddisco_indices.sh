#!/bin/bash
# filepath: test_index_changes.sh

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test configuration
TEST_DIR="./test_index_changes"
REPO_DIR=$(pwd)
BINARY_WITH_FOLDCOMP="./target/release/folddisco-with-foldcomp"
BINARY_WITHOUT_FOLDCOMP="./target/release/folddisco"
TEST_PDB_DIR="data/serine_peptidases"
TEST_FOLDCOMP_DB="data/foldcomp/example_db"
QUERY_PDB="query/4CHA.pdb"
QUERY_RESIDUES="B57,B102,C195"
NUM_THREADS=8

# Test result tracking
TEST_RESULTS_FILE="test_results.tmp"
TEST_COUNT=0

# Initialize results file
init_test_results() {
    > "$TEST_RESULTS_FILE"
}

# Helper functions
log_info() {
    echo "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo "${RED}[ERROR]${NC} $1"
}

# Test result tracking functions
record_test_result() {
    local test_name="$1"
    local status="$2"  # "PASS" or "FAIL"
    echo "$test_name:$status" >> "$TEST_RESULTS_FILE"
    TEST_COUNT=$((TEST_COUNT + 1))
}

print_results_table() {
    log_info "=== DETAILED TEST RESULTS TABLE ==="
    
    # Print table header
    printf "%-50s | %-6s\n" "Test Name" "Status"
    printf "%-50s-+-%-6s\n" "$(printf '%*s' 50 | tr ' ' '-')" "$(printf '%*s' 6 | tr ' ' '-')"
    
    local pass_count=0
    local fail_count=0
    local total_count=0
    
    # Sort test names and print results
    if [ -f "$TEST_RESULTS_FILE" ]; then
        sort "$TEST_RESULTS_FILE" | while IFS=: read -r test_name status; do
            if [ "$status" = "PASS" ]; then
                printf "%-50s | ${GREEN}%-6s${NC}\n" "$test_name" "$status"
            else
                printf "%-50s | ${RED}%-6s${NC}\n" "$test_name" "$status"
            fi
        done
        
        # Count results for summary
        pass_count=$(grep ":PASS$" "$TEST_RESULTS_FILE" | wc -l)
        fail_count=$(grep ":FAIL$" "$TEST_RESULTS_FILE" | wc -l)
        total_count=$((pass_count + fail_count))
        
        # Print summary
        printf "%-50s-+-%-6s\n" "$(printf '%*s' 50 | tr ' ' '-')" "$(printf '%*s' 6 | tr ' ' '-')"
        printf "%-50s | ${GREEN}%-6d${NC}\n" "TOTAL PASSED" "$pass_count"
        printf "%-50s | ${RED}%-6d${NC}\n" "TOTAL FAILED" "$fail_count"
        printf "%-50s | %-6d\n" "TOTAL TESTS" "$total_count"
    else
        log_warning "No test results found"
    fi
}

cleanup() {
    log_info "Cleaning up test directory..."
    if [[ -d "$TEST_DIR" ]]; then
        rm -rf "$TEST_DIR"
    fi
    # Clean up test results file
    if [[ -f "$TEST_RESULTS_FILE" ]]; then
        rm -f "$TEST_RESULTS_FILE"
    fi
}

setup_test_env() {
    log_info "Setting up test environment..."
    
    # Create test directory
    mkdir -p "$TEST_DIR"
    cd "$TEST_DIR"
    
    # Create test data if it doesn't exist
    if [[ ! -d "../$TEST_PDB_DIR" ]]; then
        log_warning "Test PDB directory not found, creating minimal test data..."
        mkdir -p "../$TEST_PDB_DIR"
        # Add code to create minimal test PDB files here
    fi
    
    # Create query PDB if it doesn't exist
    if [[ ! -f "../$QUERY_PDB" ]]; then
        log_warning "Query file not found: ../$QUERY_PDB"
        log_info "Looking for any PDB file in test directory..."
        local first_pdb=$(find "../$TEST_PDB_DIR" -name "*.pdb" | head -1)
        if [[ -n "$first_pdb" ]]; then
            QUERY_PDB="$first_pdb"
            log_info "Using: $QUERY_PDB"
        else
            log_error "No PDB files found for testing"
            exit 1
        fi
    fi
}

build_binaries() {
    log_info "Building binaries..."
    
    cd "$REPO_DIR"
    
    # Build with foldcomp feature
    log_info "Building binary with foldcomp feature..."
    cargo build --release --features foldcomp
    cp target/release/folddisco "$BINARY_WITH_FOLDCOMP" 2>/dev/null || cp target/release/folddisco "$BINARY_WITH_FOLDCOMP"
    
    # Build without foldcomp feature
    log_info "Building binary without foldcomp feature..."
    cargo build --release
    # cp target/release/folddisco "$BINARY_WITHOUT_FOLDCOMP" 2>/dev/null || cp target/release/folddisco "$BINARY_WITHOUT_FOLDCOMP"
    
    cd "$TEST_DIR"
}

test_indexing() {
    local binary="$1"
    local input_type="$2"  # "directory" or "foldcomp"
    local mode="$3"        # "id" or "big"
    local test_name="$4"
    
    log_info "Testing indexing: $test_name"
    
    local input_path
    local index_name
    
    if [[ "$input_type" == "directory" ]]; then
        input_path="../$TEST_PDB_DIR"
        index_name="idx_dir_${mode}"
    else
        input_path="../$TEST_FOLDCOMP_DB"
        index_name="idx_foldcomp_${mode}"
        
        # Skip foldcomp tests for binary without foldcomp
        if [[ "$binary" == "$BINARY_WITHOUT_FOLDCOMP" ]]; then
            log_warning "Skipping foldcomp test for binary without foldcomp feature"
            return 0
        fi
        
        # Check if foldcomp DB exists
        if [[ ! -f "$input_path" ]]; then
            log_warning "Foldcomp database not found: $input_path, skipping test"
            return 0
        fi
    fi
    
    # Test indexing
    log_info "Running: ../$binary index -p $input_path -i $index_name -m $mode -t $NUM_THREADS"

    if ../"$binary" index -p "$input_path" -i "$index_name" -m "$mode" -t "$NUM_THREADS"; then
        log_success "Indexing completed: $test_name"
        record_test_result "INDEX_$test_name" "PASS"
        
        # Verify index files exist
        local expected_files=("${index_name}.offset" "${index_name}.lookup" "${index_name}.type")
        
        # Check for value file (new format without extension)
        if [[ -f "$index_name" ]]; then
            expected_files+=("$index_name")
        elif [[ -f "${index_name}.value" ]]; then
            expected_files+=("${index_name}.value")
        else
            log_error "No value file found for $test_name"
            return 1
        fi
        
        for file in "${expected_files[@]}"; do
            if [[ ! -f "$file" ]]; then
                log_error "Missing index file: $file"
                return 1
            fi
        done
        
        # Store index info for querying tests
        echo "$index_name" >> "index_list.txt"
        
    else
        log_error "Indexing failed: $test_name"
        record_test_result "INDEX_$test_name" "FAIL"
        return 1
    fi
}

test_querying() {
    local binary="$1"
    local index_name="$2"
    local execution_location="$3"  # "repo" or "random"
    local test_name="$4"
    
    log_info "Testing querying: $test_name"
    
    local query_cmd
    local current_dir=$(pwd)
    
    if [[ "$execution_location" == "repo" ]]; then
        # Execute from repository directory
        cd "$REPO_DIR"
        query_cmd="$binary query -p $QUERY_PDB -q $QUERY_RESIDUES -i $TEST_DIR/$index_name -t $NUM_THREADS"
    else
        # Execute from random location (test directory)
        query_cmd="../$binary query -p ../$QUERY_PDB -q $QUERY_RESIDUES -i $index_name -t $NUM_THREADS"
    fi
    
    log_info "Running: $query_cmd"
    
    if $query_cmd >/dev/null 2>&1; then
        log_success "Querying completed: $test_name"
        record_test_result "QUERY_$test_name" "PASS"
    else
        log_warning "Querying failed: $test_name (continuing with other tests)"
        record_test_result "QUERY_$test_name" "FAIL"
        cd "$current_dir"
        return 0  # Return 0 to continue with other tests
    fi
    
    cd "$current_dir"
}

create_old_format_index() {
    log_info "Creating old format index for backward compatibility test..."
    
    # Create an index and then convert to old format
    local old_index="idx_old_format"

    if ../"$BINARY_WITH_FOLDCOMP" index -p "../$TEST_PDB_DIR" -i "$old_index" -m "id" -t "$NUM_THREADS"; then
        # Convert new format to old format (add .value extension)
        if [[ -f "$old_index" ]]; then
            mv "$old_index" "${old_index}.value"
            log_success "Created old format index with .value extension"
            record_test_result "INDEX_old_format_compatibility" "PASS"
            echo "$old_index" >> "index_list.txt"
            
            # Test querying with old format index immediately
            test_old_format_querying "$old_index"
        else
            log_error "Failed to create old format index"
            record_test_result "INDEX_old_format_compatibility" "FAIL"
            return 1
        fi
    else
        log_error "Failed to create base index for old format test"
        record_test_result "INDEX_old_format_compatibility" "FAIL"
        return 1
    fi
}

test_old_format_querying() {
    local old_index="$1"
    
    log_info "Testing querying with old format index (.value extension)..."
    
    # Test with foldcomp-enabled binary
    log_info "Testing old format query with foldcomp-enabled binary"
    if test_querying "$BINARY_WITH_FOLDCOMP" "$old_index" "random" "old_format_foldcomp_binary"; then
        record_test_result "QUERY_old_format_foldcomp_binary" "PASS"
    else
        record_test_result "QUERY_old_format_foldcomp_binary" "FAIL"
    fi
    
    # Test with non-foldcomp binary
    log_info "Testing old format query with non-foldcomp binary"
    if test_querying "$BINARY_WITHOUT_FOLDCOMP" "$old_index" "random" "old_format_no_foldcomp_binary"; then
        record_test_result "QUERY_old_format_no_foldcomp_binary" "PASS"
    else
        record_test_result "QUERY_old_format_no_foldcomp_binary" "FAIL"
    fi
    
    # Test from repo directory as well
    log_info "Testing old format query from repo directory"
    if test_querying "$BINARY_WITH_FOLDCOMP" "$old_index" "repo" "old_format_repo_location"; then
        record_test_result "QUERY_old_format_repo_location" "PASS"
    else
        record_test_result "QUERY_old_format_repo_location" "FAIL"
    fi
}

run_all_tests() {
    log_info "Starting comprehensive index change tests..."
    
    # Test indexing
    log_info "=== INDEXING TESTS ==="
    
    # Test with foldcomp-enabled binary
    test_indexing "$BINARY_WITH_FOLDCOMP" "directory" "id" "foldcomp_binary_dir_id"
    test_indexing "$BINARY_WITH_FOLDCOMP" "directory" "big" "foldcomp_binary_dir_big"
    test_indexing "$BINARY_WITH_FOLDCOMP" "foldcomp" "id" "foldcomp_binary_foldcomp_id"
    test_indexing "$BINARY_WITH_FOLDCOMP" "foldcomp" "big" "foldcomp_binary_foldcomp_big"
    
    # Test with non-foldcomp binary
    test_indexing "$BINARY_WITHOUT_FOLDCOMP" "directory" "id" "no_foldcomp_binary_dir_id"
    test_indexing "$BINARY_WITHOUT_FOLDCOMP" "directory" "big" "no_foldcomp_binary_dir_big"
    
    # Create old format index for backward compatibility
    create_old_format_index
    
    log_info "=== QUERYING TESTS ==="
    
    # Test querying with all created indices
    if [[ -f "index_list.txt" ]]; then
        while IFS= read -r index_name; do
            # Test from repo directory
            test_querying "$BINARY_WITH_FOLDCOMP" "$index_name" "repo" "repo_${index_name}"
            
            # Test from random location
            test_querying "$BINARY_WITH_FOLDCOMP" "$index_name" "random" "random_${index_name}"
            
            # Test with non-foldcomp binary (skip if index was created with foldcomp)
            if [[ "$index_name" != *"foldcomp"* ]]; then
                test_querying "$BINARY_WITHOUT_FOLDCOMP" "$index_name" "random" "no_foldcomp_${index_name}"
            fi
        done < "index_list.txt"
    else
        log_error "No indices were created successfully"
        return 1
    fi
}

test_default_index_path() {
    local binary="$1"
    local input_type="$2"  # "directory" or "foldcomp"
    local mode="$3"        # "id" or "big"
    local test_name="$4"
    
    log_info "Testing default index path: $test_name"
    
    local input_path
    local expected_index_name
    
    if [[ "$input_type" == "directory" ]]; then
        input_path="../$TEST_PDB_DIR"
        # Default should be input_path + "_folddisco"
        expected_index_name=$input_path"_folddisco"
    else
        input_path="../$TEST_FOLDCOMP_DB"
        # For foldcomp, default should be db_name + "_folddisco"
        expected_index_name=$input_path"_folddisco"
        
        # Skip foldcomp tests for binary without foldcomp
        if [[ "$binary" == "$BINARY_WITHOUT_FOLDCOMP" ]]; then
            log_warning "Skipping foldcomp test for binary without foldcomp feature"
            return 0
        fi
        
        # Check if foldcomp DB exists
        if [[ ! -f "$input_path" ]]; then
            log_warning "Foldcomp database not found: $input_path, skipping test"
            return 0
        fi
    fi
    
    # Test indexing WITHOUT -i flag (should use default path with _folddisco suffix)
    log_info "Running: ../$binary index -p $input_path -m $mode -t 2 (no -i flag, should create ${expected_index_name})"
    
    if ../"$binary" index -p "$input_path" -m "$mode" -t 2; then
        log_success "Default indexing completed: $test_name"
        record_test_result "DEFAULT_INDEX_$test_name" "PASS"
        
        # Verify that default indexing created files with _folddisco suffix automatically
        log_info "Verifying default index created with expected name: $expected_index_name"
        
        # Check if index files were created with expected default name
        local expected_files=("${expected_index_name}.offset" "${expected_index_name}.lookup" "${expected_index_name}.type")
        
        # Check for value file (new format without extension or old format with .value)
        if [[ -f "$expected_index_name" ]]; then
            expected_files+=("$expected_index_name")
        elif [[ -f "${expected_index_name}.value" ]]; then
            expected_files+=("${expected_index_name}.value")
        else
            log_error "No value file found for default index: $expected_index_name"
            return 1
        fi
        
        for file in "${expected_files[@]}"; do
            if [[ ! -f "$file" ]]; then
                log_error "Missing default index file: $file"
                return 1
            fi
        done
        
        # Store index info for querying tests
        echo "$expected_index_name" >> "default_index_list.txt"
        
        # Test querying with the default index
        if ! test_querying "$binary" "$expected_index_name" "random" "default_${test_name}"; then
            log_warning "Default index query test failed, but continuing..."
        fi
        
    else
        log_error "Default indexing failed: $test_name"
        record_test_result "DEFAULT_INDEX_$test_name" "FAIL"
        return 1
    fi
}

# Update the run_all_tests function to include default index path tests
run_all_tests() {
    log_info "Starting comprehensive index change tests..."
    
    # Test indexing
    log_info "=== INDEXING TESTS ==="
    
    # Test with foldcomp-enabled binary
    test_indexing "$BINARY_WITH_FOLDCOMP" "directory" "id" "foldcomp_binary_dir_id"
    test_indexing "$BINARY_WITH_FOLDCOMP" "directory" "big" "foldcomp_binary_dir_big"
    test_indexing "$BINARY_WITH_FOLDCOMP" "foldcomp" "id" "foldcomp_binary_foldcomp_id"
    test_indexing "$BINARY_WITH_FOLDCOMP" "foldcomp" "big" "foldcomp_binary_foldcomp_big"
    
    # Test with non-foldcomp binary
    test_indexing "$BINARY_WITHOUT_FOLDCOMP" "directory" "id" "no_foldcomp_binary_dir_id"
    test_indexing "$BINARY_WITHOUT_FOLDCOMP" "directory" "big" "no_foldcomp_binary_dir_big"
    
    # Create old format index for backward compatibility
    create_old_format_index
    
    log_info "=== DEFAULT INDEX PATH TESTS ==="
    
    # Test default index path behavior (no -i flag)
    test_default_index_path "$BINARY_WITH_FOLDCOMP" "directory" "id" "foldcomp_binary_dir_id_default"
    test_default_index_path "$BINARY_WITH_FOLDCOMP" "directory" "big" "foldcomp_binary_dir_big_default"
    test_default_index_path "$BINARY_WITH_FOLDCOMP" "foldcomp" "id" "foldcomp_binary_foldcomp_id_default"
    test_default_index_path "$BINARY_WITHOUT_FOLDCOMP" "directory" "id" "no_foldcomp_binary_dir_id_default"
    test_default_index_path "$BINARY_WITHOUT_FOLDCOMP" "directory" "big" "no_foldcomp_binary_dir_big_default"
    
    log_info "=== QUERYING TESTS ==="
    
    # Test querying with all created indices
    for index_list_file in "index_list.txt" "default_index_list.txt"; do
        if [[ -f "$index_list_file" ]]; then
            while IFS= read -r index_name; do
                if [[ -n "$index_name" ]]; then
                    # Test from repo directory
                    if ! test_querying "$BINARY_WITH_FOLDCOMP" "$index_name" "repo" "repo_${index_name}"; then
                        log_warning "Repo query test failed for: $index_name, but continuing..."
                    fi
                    
                    # Test from random location
                    if ! test_querying "$BINARY_WITH_FOLDCOMP" "$index_name" "random" "random_${index_name}"; then
                        log_warning "Random location query test failed for: $index_name, but continuing..."
                    fi
                    
                    # Test with non-foldcomp binary (skip if index was created with foldcomp or binary can't handle it)
                    if [[ "$index_name" != *"foldcomp"* ]]; then
                        if ! test_querying "$BINARY_WITHOUT_FOLDCOMP" "$index_name" "random" "no_foldcomp_${index_name}"; then
                            log_warning "Non-foldcomp binary query test failed for: $index_name, but continuing..."
                        fi
                    else
                        log_warning "Skipping non-foldcomp binary test for foldcomp index: $index_name"
                    fi
                fi
            done < "$index_list_file"
        fi
    done
}

# Update the main function to include shared prefix tests
main() {
    log_info "Starting comprehensive test suite for index changes..."
    
    # Initialize test results tracking
    init_test_results
    
    # Cleanup any previous test runs
    cleanup
    
    # Setup test environment
    setup_test_env
    
    # Build required binaries
    build_binaries
    
    # Run all tests
    run_all_tests
    
    # Summary
    log_info "=== TEST SUMMARY ==="
    
    # Print detailed results table
    print_results_table
    
    # Legacy summary for backward compatibility
    local failed_tests=0
    local passed_tests=0
    local total_tests=0
    
    # Count pass/fail from our tracking file
    if [ -f "$TEST_RESULTS_FILE" ]; then
        passed_tests=$(grep ":PASS$" "$TEST_RESULTS_FILE" | wc -l)
        failed_tests=$(grep ":FAIL$" "$TEST_RESULTS_FILE" | wc -l)
        total_tests=$((passed_tests + failed_tests))
    fi
    
    echo ""
    log_info "Legacy summary:"
    log_info "Total tests executed: $total_tests"
    log_success "Tests passed: $passed_tests"
    
    if [[ $failed_tests -gt 0 ]]; then
        log_error "Tests failed: $failed_tests"
        log_error "Check the detailed results table above for failure details"
        exit 1
    else
        log_success "All tests passed successfully!"
    fi
    
    # Cleanup
    cd "$REPO_DIR"
    log_info "Test completed. Test files are in $TEST_DIR/"
    log_info "Run 'rm -rf $TEST_DIR' to clean up test files"
}

# Handle script interruption
trap cleanup EXIT

# Run main function
main "$@"