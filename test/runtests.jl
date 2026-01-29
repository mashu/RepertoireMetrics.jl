using RepertoireMetrics
using Test
using DataFrames

@testset "RepertoireMetrics.jl" begin
    
    @testset "Types - Lineage Definitions" begin
        # LineageIDDefinition
        lid = LineageIDDefinition()
        @test lid.column == :lineage_id
        
        lid2 = LineageIDDefinition(:clone_id)
        @test lid2.column == :clone_id
        
        # VJCdr3Definition
        vjc = VJCdr3Definition()
        @test vjc.v_column == :v_call
        @test vjc.j_column == :j_call
        @test vjc.cdr3_column == :cdr3
        @test vjc.use_first_allele == true
        
        vjc2 = VJCdr3Definition(use_first_allele=false)
        @test vjc2.use_first_allele == false
        
        # CustomDefinition
        custom = CustomDefinition(row -> row.v_call)
        @test custom.key_func isa Function
    end
    
    @testset "first_allele" begin
        @test first_allele("IGHV1-2*01") == "IGHV1-2*01"
        @test first_allele("IGHV1-2*01,IGHV1-2*02") == "IGHV1-2*01"
        @test first_allele("IGHV1-2*01, IGHV1-2*02") == "IGHV1-2*01"
        @test first_allele("") == ""
        @test first_allele(missing) === missing
    end
    
    @testset "Repertoire Construction" begin
        # Basic construction
        rep = Repertoire([100, 50, 25, 10, 5])
        @test richness(rep) == 5
        @test total_count(rep) == 190
        @test length(rep) == 5
        @test !isempty(rep)
        
        # Counts should be sorted descending
        @test counts(rep) == [100, 50, 25, 10, 5]
        
        # With custom IDs
        rep2 = Repertoire(
            [50, 100, 25],  # Will be sorted
            ["clone_b", "clone_a", "clone_c"];
            donor_id = "Donor1"
        )
        @test counts(rep2) == [100, 50, 25]
        @test lineage_ids(rep2) == ["clone_a", "clone_b", "clone_c"]
        @test donor_id(rep2) == "Donor1"
        
        # Empty repertoire
        empty_rep = Repertoire(Int[])
        @test richness(empty_rep) == 0
        @test isempty(empty_rep)
        @test total_count(empty_rep) == 0
    end
    
    @testset "Frequencies" begin
        rep = Repertoire([50, 30, 20])
        freqs = frequencies(rep)
        @test length(freqs) == 3
        @test freqs ≈ [0.5, 0.3, 0.2]
        @test sum(freqs) ≈ 1.0
        
        # Empty repertoire
        empty_rep = Repertoire(Int[])
        @test frequencies(empty_rep) == Float64[]
    end
    
    @testset "Shannon Entropy" begin
        # Perfect evenness with 4 species: H = log(4)
        rep_even = Repertoire([25, 25, 25, 25])
        @test shannon_entropy(rep_even) ≈ log(4) atol=1e-10
        
        # Single species: H = 0
        rep_single = Repertoire([100])
        @test shannon_entropy(rep_single) ≈ 0.0 atol=1e-10
        
        # Two species, 50-50: H = log(2)
        rep_two = Repertoire([50, 50])
        @test shannon_entropy(rep_two) ≈ log(2) atol=1e-10
        
        # Known calculation: p = [0.5, 0.3, 0.2]
        rep_known = Repertoire([50, 30, 20])
        expected_H = -0.5*log(0.5) - 0.3*log(0.3) - 0.2*log(0.2)
        @test shannon_entropy(rep_known) ≈ expected_H atol=1e-10
    end
    
    @testset "Simpson Index" begin
        # Perfect evenness with n species: D = 1/n
        rep_even = Repertoire([25, 25, 25, 25])
        @test simpson_index(rep_even) ≈ 0.25 atol=1e-10
        
        # Single species: D = 1
        rep_single = Repertoire([100])
        @test simpson_index(rep_single) ≈ 1.0 atol=1e-10
        
        # Known calculation: D = Σp²
        rep_known = Repertoire([50, 30, 20])
        expected_D = 0.5^2 + 0.3^2 + 0.2^2
        @test simpson_index(rep_known) ≈ expected_D atol=1e-10
        
        # Simpson diversity = 1 - D
        @test simpson_diversity(rep_known) ≈ (1 - expected_D) atol=1e-10
        
        # Inverse Simpson
        @test inverse_simpson(rep_known) ≈ (1 / expected_D) atol=1e-10
    end
    
    @testset "Clonality and Evenness" begin
        # Perfect evenness: clonality = 0, evenness = 1
        rep_even = Repertoire([25, 25, 25, 25])
        @test clonality(rep_even) ≈ 0.0 atol=1e-10
        @test evenness(rep_even) ≈ 1.0 atol=1e-10
        
        # Single species: clonality = 1
        rep_single = Repertoire([100])
        @test clonality(rep_single) ≈ 1.0 atol=1e-10
    end
    
    @testset "Berger-Parker Index" begin
        rep = Repertoire([50, 30, 20])
        @test berger_parker_index(rep) ≈ 0.5 atol=1e-10
        
        rep2 = Repertoire([80, 10, 10])
        @test berger_parker_index(rep2) ≈ 0.8 atol=1e-10
    end
    
    @testset "Gini Coefficient" begin
        # Perfect equality: Gini = 0
        rep_equal = Repertoire([25, 25, 25, 25])
        @test gini_coefficient(rep_equal) ≈ 0.0 atol=1e-10
        
        # High inequality (for 2 species, max Gini is 0.5)
        rep_unequal = Repertoire([99, 1])
        @test gini_coefficient(rep_unequal) > 0.4
        
        # More species, higher inequality possible
        rep_very_unequal = Repertoire([1000, 1, 1, 1, 1])
        @test gini_coefficient(rep_very_unequal) > 0.7
        
        # Single species
        rep_single = Repertoire([100])
        @test gini_coefficient(rep_single) ≈ 0.0 atol=1e-10
    end
    
    @testset "D50" begin
        # 50% threshold
        rep = Repertoire([50, 30, 10, 5, 5])  # Total = 100
        @test d50(rep) == 1  # First clone = 50% exactly
        
        rep2 = Repertoire([40, 30, 20, 10])  # Total = 100
        @test d50(rep2) == 2  # 40 + 30 = 70 >= 50
        
        rep3 = Repertoire([20, 20, 20, 20, 20])  # Even distribution
        @test d50(rep3) == 3  # 20+20+20 = 60 >= 50
    end
    
    @testset "Chao1" begin
        # No singletons or doubletons: Chao1 = observed richness
        rep_no_rare = Repertoire([100, 50, 30])
        @test chao1(rep_no_rare) ≈ 3.0 atol=1e-10
        
        # With singletons and doubletons
        rep_with_rare = Repertoire([100, 50, 2, 2, 1, 1, 1])  # f1=3, f2=2
        expected_chao1 = 7 + (3*3)/(2*2)
        @test chao1(rep_with_rare) ≈ expected_chao1 atol=1e-10
    end
    
    @testset "Hill Numbers" begin
        rep = Repertoire([50, 30, 20])
        freqs = frequencies(rep)
        
        # q=0: Richness
        h0 = hill_number(rep, 0)
        @test h0.value ≈ 3.0 atol=1e-10
        
        # q=1: exp(Shannon entropy)
        h1 = hill_number(rep, 1)
        @test h1.value ≈ exp(shannon_entropy(rep)) atol=1e-10
        
        # q=2: Inverse Simpson
        h2 = hill_number(rep, 2)
        @test h2.value ≈ inverse_simpson(rep) atol=1e-10
    end
    
    @testset "compute_metrics - default all" begin
        rep = Repertoire([50, 30, 20])
        metrics = compute_metrics(rep)
        
        @test metrics isa Metrics
        @test metrics.richness == 3
        @test metrics.total_count == 100
        @test metrics.shannon_entropy ≈ shannon_entropy(rep) atol=1e-10
        @test metrics.simpson_index ≈ simpson_index(rep) atol=1e-10
        @test metrics.clonality ≈ clonality(rep) atol=1e-10
        @test metrics.d50 == d50(rep)
        @test metrics.chao1 ≈ chao1(rep) atol=1e-10
    end
    
    @testset "DataFrame Conversion" begin
        # Create test DataFrame with VJCdr3
        df = DataFrame(
            v_call = ["IGHV1-2*01", "IGHV1-2*01", "IGHV3-21*01", "IGHV3-21*01"],
            j_call = ["IGHJ4*02", "IGHJ4*02", "IGHJ6*01", "IGHJ6*01"],
            cdr3 = ["CARDYW", "CARDYW", "CAKGW", "CARDYW"],
            count = [10, 5, 20, 3]
        )
        
        # Using VJCdr3Definition should aggregate by V-J-CDR3
        rep = repertoire_from_dataframe(df, VJCdr3Definition())
        @test richness(rep) == 3  # Three unique V-J-CDR3 combinations
        @test total_count(rep) == 38
        
        # Test without count column (each row = 1)
        rep_no_count = repertoire_from_dataframe(
            df, VJCdr3Definition();
            count_column = nothing
        )
        @test total_count(rep_no_count) == 4
        
        # Using LineageIDDefinition
        df_with_id = DataFrame(
            lineage_id = ["L1", "L1", "L2", "L3"],
            count = [10, 5, 20, 3]
        )
        rep_id = repertoire_from_dataframe(df_with_id, LineageIDDefinition())
        @test richness(rep_id) == 3
        @test total_count(rep_id) == 38
    end
    
    @testset "RepertoireCollection" begin
        rep1 = Repertoire([100, 50], donor_id="D1")
        rep2 = Repertoire([80, 60, 40], donor_id="D2")
        
        collection = RepertoireCollection([rep1, rep2])
        @test length(collection) == 2
        @test collection.donor_ids == ["D1", "D2"]
        
        # Iteration
        donors_found = String[]
        for rep in collection
            push!(donors_found, donor_id(rep))
        end
        @test donors_found == ["D1", "D2"]
        
        # Indexing
        @test donor_id(collection[1]) == "D1"
        @test donor_id(collection[2]) == "D2"
        
        # Compute metrics for collection
        all_metrics = compute_metrics(collection)
        @test length(all_metrics) == 2
        @test all_metrics[1].richness == 2
        @test all_metrics[2].richness == 3
    end
    
    @testset "metrics_to_dataframe" begin
        rep = Repertoire([50, 30, 20], donor_id="TestDonor")
        metrics = compute_metrics(rep)
        
        df = metrics_to_dataframe(metrics, "TestDonor")
        @test nrow(df) == 1
        @test df.donor_id[1] == "TestDonor"
        @test df.richness[1] == 3
        @test df.clonality[1] ≈ clonality(rep) atol=1e-10
    end
    
    @testset "Type Stability" begin
        # Test that key functions return consistent types
        rep_int = Repertoire([100, 50, 25])
        rep_float = Repertoire([100.0, 50.0, 25.0])
        
        @test typeof(richness(rep_int)) == Int
        @test typeof(richness(rep_float)) == Int
        
        @test typeof(shannon_entropy(rep_int)) == Float64
        @test typeof(shannon_entropy(rep_float)) == Float64
        
        @test typeof(frequencies(rep_int)) == Vector{Float64}
        @test typeof(frequencies(rep_float)) == Vector{Float64}
        
        @test typeof(compute_metrics(rep_int)) == Metrics
        @test typeof(compute_metrics(rep_float)) == Metrics
    end
    
    @testset "Edge Cases" begin
        # Empty repertoire
        empty_rep = Repertoire(Int[])
        metrics = compute_metrics(empty_rep)
        @test metrics.richness == 0
        @test metrics.total_count == 0
        @test metrics.shannon_entropy == 0.0
        
        # Single element
        single_rep = Repertoire([100])
        metrics = compute_metrics(single_rep)
        @test metrics.richness == 1
        @test metrics.shannon_entropy == 0.0
        @test metrics.clonality == 1.0
        @test metrics.d50 == 1
    end
    
    @testset "Lineage Key Extraction" begin
        # Create a row-like named tuple
        row = (v_call = "IGHV1-2*01,IGHV1-2*02", j_call = "IGHJ4*02", cdr3 = "CARDYW", lineage_id = "L1")
        
        # VJCdr3 with first allele
        strategy1 = VJCdr3Definition()
        key1 = lineage_key(strategy1, row)
        @test key1 == ("IGHV1-2*01", "IGHJ4*02", "CARDYW")
        
        # VJCdr3 without first allele
        strategy2 = VJCdr3Definition(use_first_allele=false)
        key2 = lineage_key(strategy2, row)
        @test key2 == ("IGHV1-2*01,IGHV1-2*02", "IGHJ4*02", "CARDYW")
        
        # LineageID
        strategy3 = LineageIDDefinition()
        key3 = lineage_key(strategy3, row)
        @test key3 == "L1"
        
        # Custom
        strategy4 = CustomDefinition(r -> r.cdr3)
        key4 = lineage_key(strategy4, row)
        @test key4 == "CARDYW"
    end
    
    @testset "Composable MetricSet" begin
        rep = Repertoire([50, 30, 20])
        
        # Test + operator combining metrics
        metrics_set = ShannonEntropy() + Clonality()
        @test metrics_set isa MetricSet
        @test length(metrics_set) == 2
        
        # Test computing selected metrics
        result = compute_metrics(rep, metrics_set)
        @test result isa Metrics
        @test result.shannon_entropy ≈ shannon_entropy(rep) atol=1e-10
        @test result.clonality ≈ clonality(rep) atol=1e-10
        
        # Missing metrics return missing
        @test result.d50 === missing
        
        # Test combining more metrics
        extended = metrics_set + D50() + Richness()
        @test length(extended) == 4
        
        result2 = compute_metrics(rep, extended)
        @test result2.d50 == d50(rep)
        @test result2.richness == richness(rep)
        
        # Test predefined metric sets
        div_result = compute_metrics(rep, DIVERSITY_METRICS)
        @test div_result.shannon_entropy ≈ shannon_entropy(rep) atol=1e-10
        @test div_result.evenness ≈ evenness(rep) atol=1e-10
        
        clon_result = compute_metrics(rep, CLONALITY_METRICS)
        @test clon_result.clonality ≈ clonality(rep) atol=1e-10
        @test clon_result.gini_coefficient ≈ gini_coefficient(rep) atol=1e-10
        
        # Test single metric
        single = compute_metrics(rep, ShannonEntropy())
        @test single.shannon_entropy ≈ shannon_entropy(rep) atol=1e-10
        
        # Test compute_metric for individual metrics
        @test compute_metric(rep, ShannonEntropy()) ≈ shannon_entropy(rep) atol=1e-10
        @test compute_metric(rep, Clonality()) ≈ clonality(rep) atol=1e-10
        @test compute_metric(rep, Richness()) == richness(rep)
    end
    
    @testset "Depth metric and ROBUST_METRICS" begin
        rep = Repertoire([100, 50, 25, 10, 5])
        
        # Test Depth metric
        metrics = compute_metrics(rep, Depth())
        @test metrics.depth == 190.0  # Float64 conversion
        
        # Test Depth equals TotalCount
        metrics2 = compute_metrics(rep, Depth() + TotalCount())
        @test metrics2.depth == metrics2.total_count
        
        # Test ROBUST_METRICS includes depth
        robust = compute_metrics(rep, ROBUST_METRICS)
        @test haskey(robust.values, :depth)
        @test haskey(robust.values, :simpson_diversity)
        @test haskey(robust.values, :clonality)
        @test haskey(robust.values, :gini_coefficient)
    end
    
    @testset "MetricSet with Collection" begin
        rep1 = Repertoire([100, 50], donor_id="D1")
        rep2 = Repertoire([80, 60, 40], donor_id="D2")
        collection = RepertoireCollection([rep1, rep2])
        
        # Compute selected metrics for collection
        metrics_set = ShannonEntropy() + Clonality()
        results = compute_metrics(collection, metrics_set)
        
        @test length(results) == 2
        @test results[1].shannon_entropy ≈ shannon_entropy(rep1) atol=1e-10
        @test results[2].clonality ≈ clonality(rep2) atol=1e-10
        
        # Convert to DataFrame
        df = metrics_to_dataframe(collection, results)
        @test nrow(df) == 2
        @test :shannon_entropy in propertynames(df)
        @test :clonality in propertynames(df)
    end
    
    @testset "Length Statistics (Composable)" begin
        # Create test DataFrame with sequences (nucleotide, length divisible by 3)
        df = DataFrame(
            v_call = ["IGHV1-1*01", "IGHV1-2*01", "IGHV1-3*01", "IGHV1-4*01"],
            j_call = ["IGHJ1*01", "IGHJ1*01", "IGHJ2*01", "IGHJ2*01"],
            cdr3 = ["TGTGCCAGG", "TGTGCCAGGAAA", "TGTGCC", "TGTGCCAGGAAATTT"],  # lengths 9,12,6,15 nt -> 3,4,2,5 aa
            count = [10, 5, 3, 2]
        )
        
        # Test compute_length_stats directly (unweighted)
        stats = compute_length_stats(df; length_column=:cdr3)
        @test stats.n_sequences == 4
        @test stats.min_length == 2
        @test stats.max_length == 5
        @test stats.mean_length ≈ (3 + 4 + 2 + 5) / 4 atol=1e-10
        @test stats.column == :cdr3
        
        # Weighted statistics
        stats_w = compute_length_stats(df; length_column=:cdr3, count_column=:count)
        @test stats_w.n_sequences == 20  # 10+5+3+2
        # Weighted mean: (3*10 + 4*5 + 2*3 + 5*2) / 20 = (30+20+6+10)/20 = 66/20 = 3.3
        @test stats_w.mean_length ≈ 3.3 atol=1e-10
        
        # Length distribution
        dist = length_distribution(df; length_column=:cdr3, count_column=:count)
        @test dist[3] == 10  # length 3 has count 10
        @test dist[4] == 5   # length 4 has count 5
        @test dist[2] == 3   # length 2 has count 3
        @test dist[5] == 2   # length 5 has count 2
        
        # Test with amino acid sequences
        df_aa = DataFrame(
            v_call = ["IGHV1-1*01", "IGHV1-2*01", "IGHV1-3*01", "IGHV1-4*01"],
            j_call = ["IGHJ1*01", "IGHJ1*01", "IGHJ2*01", "IGHJ2*01"],
            cdr3 = ["AAA", "AAAAAA", "AAAAAAAAA", "AAAAAAAAAAAA"],  # just for lineage
            cdr3_aa = ["CAR", "CARD", "CA", "CARDY"],  # 3, 4, 2, 5
            count = [10, 5, 3, 2]
        )
        stats_aa = compute_length_stats(df_aa; length_column=:cdr3_aa, use_aa=true)
        @test stats_aa.min_length == 2
        @test stats_aa.max_length == 5
        
        # Test composable integration - create repertoire with length stats
        rep = repertoire_from_dataframe(df, VJCdr3Definition(); 
            count_column=:count, length_column=:cdr3)
        
        @test has_length_stats(rep)
        @test mean_length(rep) ≈ 3.3 atol=1e-10
        @test min_length(rep) == 2
        @test max_length(rep) == 5
        
        # Composable metrics
        metrics = compute_metrics(rep, MeanLength() + MedianLength())
        @test metrics.mean_length ≈ 3.3 atol=1e-10
        @test haskey(metrics.values, :median_length)
        
        # All length metrics at once
        metrics_all = compute_metrics(rep, LENGTH_METRICS)
        @test haskey(metrics_all.values, :mean_length)
        @test haskey(metrics_all.values, :std_length)
        @test haskey(metrics_all.values, :min_length)
        @test haskey(metrics_all.values, :max_length)
        
        # Mix with other metrics
        mixed = compute_metrics(rep, Richness() + MeanLength() + ShannonEntropy())
        @test haskey(mixed.values, :richness)
        @test haskey(mixed.values, :mean_length)
        @test haskey(mixed.values, :shannon_entropy)
        
        # Empty DataFrame
        empty_df = DataFrame(cdr3 = String[], v_call = String[], j_call = String[])
        empty_stats = compute_length_stats(empty_df; length_column=:cdr3)
        @test empty_stats.n_sequences == 0
        
        # Repertoire without length stats should error on access
        rep_no_length = Repertoire([100, 50, 25])
        @test !has_length_stats(rep_no_length)
        @test_throws ErrorException mean_length(rep_no_length)
    end
    
end
