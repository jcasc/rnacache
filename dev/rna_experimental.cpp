using tax_distribution = std::unordered_map<const taxon*, double>;
using atomic_tax_distribution = std::unordered_map<const taxon*, std::atomic<double>>;
tax_distribution ground_distribution(const database& db, const rna_ground_truth& mappings) {
    tax_distribution rho, alpha;
    double alpha_sum = 0;

    for (const auto& map: mappings)
        rho[map]++;
   
    for (auto& i: rho)
        i.second /= mappings.size();
    
    for (const auto& it: db.target_taxa()) {
        int target_length = (it.source().windows-1)*db.target_window_stride()+db.target_window_size();
        alpha[&it] = rho[&it]/target_length;
        alpha_sum += alpha[&it];
    }

    for (auto& i: alpha)
        i.second /= alpha_sum;

    return alpha;
}

tax_distribution ml_distribution_parallel(const database& db, const query_options& opt, const query_matches& mappings, size_t iterations=10) {
    auto rho = std::make_unique<atomic_tax_distribution>();

    for (const auto& it: db.target_taxa())
        rho->emplace(&it, 1.0 / db.target_count());
    
    // Iteration
    for (size_t iter = 0; iter < iterations; iter++) {
        std::cerr << "EM Iteration: " << iter << std::endl;
        auto tax_count = std::make_unique<atomic_tax_distribution>();
        for (const auto& it: db.target_taxa())
            tax_count->emplace(&it, 0);

        std::vector<std::thread> threads;
        threads.reserve(opt.process.numThreads);
        for (size_t p = 0; p < opt.process.numThreads; ++p)
            threads.emplace_back([&,p]{
                size_t start = p*(mappings.size()/opt.process.numThreads)+std::min(p, mappings.size()%opt.process.numThreads);
                size_t end = (p+1)*(mappings.size()/opt.process.numThreads)+std::min(p+1, mappings.size()%opt.process.numThreads);
                for (size_t i = start; i<end; ++i) {
                    tax_distribution z;
                    double z_sum = 0;
                    for (const taxon* tax: mappings[i]) {
                        int target_length = (tax->source().windows-1)*db.target_window_stride()+db.target_window_size();
                        z[tax] = (*rho)[tax] / (target_length - opt.classify.insertSizeMax + 1);
                        z_sum += z[tax];
                    }
                    for (const auto& z_: z) {
                        double tmp = (*tax_count)[z_.first];
                        double summand = (z_.second/z_sum)/mappings.size();
                        while (!(*tax_count)[z_.first].compare_exchange_weak(tmp, tmp + summand));
                    }
                }
            });
        for (auto& t: threads) t.join();
        rho = std::move(tax_count);
    }

    tax_distribution alpha;
    double alpha_sum = 0;
    for (const auto& it: db.target_taxa()) {
        int target_length = (it.source().windows-1)*db.target_window_stride()+db.target_window_size();
        alpha[&it] = (*rho)[&it]/target_length;
        alpha_sum += alpha[&it];
    }

    for (auto& i: alpha)
        i.second /= alpha_sum;
    
    return alpha;
}

double spearman_rho(const database& db, const tax_distribution& x, const tax_distribution& y) {
    std::vector<const taxon*> xs, ys;
    tax_distribution xr, yr;
    for (const auto& it: db.target_taxa()) {
        xs.emplace_back(&it);
        ys.emplace_back(&it);
    }

    std::sort(xs.begin(), xs.end(), [&](const taxon* lhs, const taxon* rhs) {return x.at(lhs) < x.at(rhs);});
    std::sort(ys.begin(), ys.end(), [&](const taxon* lhs, const taxon* rhs) {return y.at(lhs) < y.at(rhs);});

    size_t sum_ranks_x = 0, sum_ranks_y = 0;
    std::set<const taxon*> x_tmp{xs.front()};
    std::set<const taxon*> y_tmp{ys.front()};
    for (size_t i = 1; i < xs.size(); i++) {
        if (x.at(xs[i]) == x.at(xs[i-1])) {
            sum_ranks_x += i;
            x_tmp.emplace(xs[i]);
        }
        else {
            for (const auto& j: x_tmp)
                xr[j] = double(sum_ranks_x) / x_tmp.size();
            sum_ranks_x = i;
            x_tmp = {xs[i]};
        }
        if (y.at(ys[i]) == y.at(ys[i-1])) {
            sum_ranks_y += i;
            y_tmp.emplace(ys[i]);
        }
        else {
            for (const auto& j: y_tmp)
                yr[j] = double(sum_ranks_y) / y_tmp.size();
            sum_ranks_y = i;
            y_tmp = {ys[i]};
        }
    }
    for (const auto& j: x_tmp)
        xr[j] = double(sum_ranks_x) / x_tmp.size();
    for (const auto& j: y_tmp)
        yr[j] = double(sum_ranks_y) / y_tmp.size();

    // for (const auto& i: xr)
    //     std::cerr << i.second << " " << yr[i.first] << std::endl;

    double mean_x = 0, mean_y = 0, sd_x = 0, sd_y = 0, cov_xy = 0;
    for (const auto& i: xr) {
        mean_x += i.second;
        mean_y += yr[i.first];
    }
    mean_x /= xr.size();
    mean_y /= yr.size();
    
    for (const auto& i: xr) {
        sd_x += (i.second - mean_x)*(i.second - mean_x);
        sd_y += (yr[i.first] - mean_y)*(yr[i.first] - mean_y);
        cov_xy += (i.second - mean_x)*(yr[i.first] - mean_y);
    }
    sd_x = std::sqrt(sd_x/(xr.size()-1));
    sd_y = std::sqrt(sd_y/(yr.size()-1));
    cov_xy /= (xr.size()-1);

    return (cov_xy/(sd_x*sd_y));
}


void prefix_sum(std::vector<size_t> dims, std::unordered_map<size_t, size_t> &data) {
	std::vector<size_t> steps(1,1);
	for (size_t i = 0; i<dims.size()-1; i++)
		steps.push_back(steps[i]*dims[i]);

	for (size_t dim = 0; dim<dims.size(); dim++) {
		//std::cout << "DIM: " << dim << std::endl;
		std::vector<size_t> pos(dims.size(), 0);
		bool b = false;
		while(!b) {
			size_t offset = 0;
			for (size_t i = 0; i<pos.size(); i++)
				offset+=pos[i]*steps[i];
			for (size_t i = dims[dim]-2; i<dims[dim]-1; i--)
				data[offset+i*steps[dim]]+=data[offset+(i+1)*steps[dim]];
				//std::cout << offset+i*steps[dim] << " += " << offset+(i+1)*steps[dim] << std::endl;
			b = true;
			for (size_t i = pos.size()-1; i<pos.size(); i--) {
				if (i==dim) continue;
				if (pos[i]<dims[i]-1) {
					pos[i]++;
					for (size_t j=i+1; j<pos.size(); j++)
						pos[j]=0;
					b=false;
					break;
				}
			}
		}
		//std::cout << std::endl;
	}
}

/*************************************************************************//**
 *
 * @brief parameter search for rna mapping
 *
 *****************************************************************************/
void map_queries_to_targets_rna_search(
    const vector<string>& infiles,
    const database& db, const query_options& opt)
{    
    struct eval_result {
    // parameters
    size_t hit_thresh = 0;
    double cutoff = 0;
    double cov_stat = 0;

    // results
    size_t reads_total = 0;

    size_t matches_total = 0; 
    size_t reads_aligned = 0;
    size_t origin_mapped = 0;
    size_t correctly_rejected = 0;
    
    };  

    query_matches mappings_;
    taxon_list truth_;
    std::vector<std::vector<uint16_t>> hitcounts_;
    matches_per_target_param pcoverage_;

    const auto makeBatchBuffer = [] { return mappings_buffer(); };

    const auto processQuery = [&] (mappings_buffer& buf,
        sequence_query&& query, match_locations& allhits)
    {
        if(query.empty()) return;
        
        buf.truth.emplace_back(ground_truth_taxon(db, query.header));
                
        taxon_list cls;
        std::vector<uint16_t> cls_hits;
        auto cands = make_classification_candidates(db, opt.classify, query, allhits);
        
        for (const auto& cand: cands) {
            cls.emplace_back(cand.tax);
            cls_hits.emplace_back(cand.hits);
        }

        buf.mappings.emplace_back(std::move(cls));
        buf.hitcounts.emplace_back(std::move(cls_hits));
        buf.pcoverage.insert(allhits, cands);
    };

    //runs before a batch buffer is discarded
    const auto finalizeBatch = [&] (mappings_buffer&& buf) {
        truth_.insert(truth_.end(), buf.truth.begin(), buf.truth.end());
        mappings_.insert(mappings_.end(),
                          std::make_move_iterator(buf.mappings.begin()),
                          std::make_move_iterator(buf.mappings.end() ));
        hitcounts_.insert(hitcounts_.end(),
                          std::make_move_iterator(buf.hitcounts.begin()),
                          std::make_move_iterator(buf.hitcounts.end() ));
        pcoverage_.merge(std::move(buf.pcoverage));
        std::cerr << mappings_.size() << std::endl;
    };

    //runs if something needs to be appended to the output
    const auto appendToOutput = [] (const std::string&) {};

    //run (parallel) database queries according to processing options
    query_database(infiles, db, opt.process,
                   makeBatchBuffer, processQuery, finalizeBatch,
                   appendToOutput);

    std::vector<std::pair<uint16_t, double>> settings;
    // for (uint8_t rel = 0; rel <=100; rel+=1)
    //         settings.emplace_back(5, double(rel)/100);
    settings.emplace_back(5, 0.8);
    
    std::vector<std::thread> threads;
    std::mutex eval_lock;
    std::atomic<int> queue{0};
    struct cmp {
        bool operator() (const eval_result& lhs, const eval_result& rhs) const noexcept{
            if (lhs.matches_total == rhs.matches_total)
                return lhs.origin_mapped > rhs.origin_mapped;
            else
                return lhs.matches_total < rhs.matches_total;
        }
    };
    std::multiset<eval_result, cmp> global_results;

    for (size_t p = 0; p < opt.process.numThreads; ++p) {
        threads.emplace_back([&,p]{
            for (size_t queue_; (queue_ = queue++) < settings.size();) {
                std::cerr << p << " claimed " << queue_ << std::endl;
                uint16_t abs = settings[queue_].first;
                double rel = settings[queue_].second;
                constexpr size_t N = 200;
                
                // ~sparse vectors~
                std::unordered_map<size_t, size_t> total_matches, origin_mapped, aligned, correctly_rejected;
                
                std::unordered_map<const taxon* , double> cov_cache;

                for (size_t qid = 0; qid < mappings_.size(); ++qid) {
                    auto cands = mappings_[qid];
                    auto cand_hits = hitcounts_[qid];
                    
                    if (cands.empty()) {
                        if (truth_[qid] == nullptr)
                            correctly_rejected[0]++;
                        continue;
                    }
                    
                    size_t max_hits = *std::max_element(cand_hits.cbegin(), cand_hits.cend());
                    
                    double max_cov = 0;
                    for (size_t cid = 0; cid < cands.size(); ++cid) {
                        if (cand_hits[cid] < abs || double(cand_hits[cid])/max_hits < rel)
                            continue;

                        if (!cov_cache.count(cands[cid]))
                            cov_cache[cands[cid]] = pcoverage_.coverage(cands[cid], abs, rel);
                        
                        max_cov = std::max(max_cov, cov_cache[cands[cid]]);
                    }
                    
                    double max_cov_stat = 0;
                    double min_cov_stat = 1.0;
                    for (size_t cid = 0; cid < cands.size(); ++cid) {
                        if (cand_hits[cid] < abs || double(cand_hits[cid])/max_hits < rel)
                            continue;
                        
                        double cov_stat = 1.0;
                        if (max_cov > 0.0)
                            cov_stat = cov_cache[cands[cid]]/max_cov;

                        max_cov_stat = std::max(max_cov_stat, cov_stat);
                        min_cov_stat = std::min(min_cov_stat, cov_stat);
                                                
                        total_matches[cov_stat*N]++;
                        if (((cov_stat*N)/N < 0.8 && cov_stat >= 0.8) ||
                           ((cov_stat*N)/N >= 0.8 && cov_stat < 0.8)) {
                            std::cout << "BRO YOU JUST POSTED CRINGE BRO" << std::endl;
                            std::cout << std::setprecision(100) << 0.8 << " " << cov_stat << " " << (cov_stat*N)/N << std::endl;
                        }

                        if (truth_[qid] == cands[cid])
                            origin_mapped[cov_stat*N]++;
                    }
                    aligned[max_cov_stat*N]++;
                    if (truth_[qid] == nullptr)
                        correctly_rejected[std::ceil(min_cov_stat*N)]++;
                }

                // prefix_sums
                for (size_t i = N-1; i < N; --i) {
                    total_matches[i] += total_matches[i+1];
                    origin_mapped[i] += origin_mapped[i+1];
                    aligned[i] += aligned[i+1];
                    correctly_rejected[N-i] += correctly_rejected[N-i-1];
                }
	        
                std::cerr << p << " finished " << queue_ << std::endl;
                std::lock_guard<std::mutex> guard(eval_lock);
                for (size_t i = 0; i <= N; ++i) {
                    // filter here if neccessary
                    global_results.insert({abs, rel, double(i)/N, mappings_.size(), total_matches[i], aligned[i], origin_mapped[i], correctly_rejected[i]});
                }
	            std::cerr << p << " finalized " << queue_ << std::endl;
            }
        });
    }

    for (std::thread& p: threads) p.join();
    std::cerr << global_results.size() << " results" << std::endl;
    std::cout << "hit_threshold,hit_cutoff,relative_coverage,reads,matches,aligned,TP,TN" << std::endl;
    for (const auto& i: global_results) {
        std::cout << i.hit_thresh << "," << i.cutoff << "," << i.cov_stat << "," 
                << i.reads_total << "," << i.matches_total << "," << i.reads_aligned << "," << i.origin_mapped << "," << i.correctly_rejected << std::endl;
    }
    std::cout << std::endl;
}

//// NEW
//// 1 PASS STUFF


//template <class Consumer>
//void do_in_parallel(size_t num_threads, size_t num_indices, Consumer&& consume) {
//    std::vector<std::thread> threads;
//    for (size_t thid = 0; thid < num_threads; ++thid)
//        threads.emplace_back([&,thid]{
//            size_t start = thid * (num_indices/num_threads) + std::min(thid, num_indices%num_threads);
//            size_t end = (thid+1) * (num_indices/num_threads) + std::min(thid+1, num_indices%num_threads);
//            for (size_t i = start; i < end; ++i)
//                consume(i);
//        });
//    for (auto& t: threads) t.join();
//
//}

// needed for 1pass
// /*************************************************************************//**
//  *
//  * @brief for each query applies coverage filter to list of taxon matches
//  *
//  *****************************************************************************/
// void coverage_filter(query_matches& matches,
//                      const matches_per_target_light& coverage,
//                      const query_options& opt)
// {
//     do_in_parallel(opt.process.numThreads, matches.size(), [&](std::size_t i) {
//         coverage_filter(matches[i], coverage, opt);
//     });
// }


// needed for 1pass
// /*************************************************************************//**
//  *
//  * @brief
//  *
//  *****************************************************************************/
// void evaluate_classification(
//     const database& db,
//     const query_options& opt,
//     const query_matches& matches,
//     const taxon_list& truth,
//     rna_mapping_statistics& statistics)
// {
//     if(opt.output.evaluate.determineGroundTruth) {
//         do_in_parallel(opt.performance.numThreads, matches.size(), [&](size_t i) {
//             evaluate_classification(opt.evaluate, truth[i], matches[i], statistics);
//         });
//     } else {
//         do_in_parallel(opt.performance.numThreads, matches.size(), [&](size_t i) {
//             evaluate_classification(opt.evaluate, nullptr, matches[i], statistics);
//         });
//     }
// }


// /*************************************************************************//**
//  *
//  * @brief transcriptome classification scheme; 1-pass variant;
//  *        needs more memory, but is potentially faster
//  *
//  *****************************************************************************/
// void map_queries_to_targets_1pass(
//     const vector<string>& infiles,
//     const database& db,
//     const query_options& opt,
//     classification_results& results)
// {
//     //global target -> (hit_windows set) map
//     matches_per_target_light coverage_;
//     query_matches mappings_;
//     taxon_list truth_;

//     //input queries are divided into batches;
//     //each batch might be processed by a different thread;
//     //the following 4 lambdas define actions that should be performed
//     //on such a batch and its associated buffer;
//     //the batch buffer can be used to cache intermediate results

//     //creates an empty batch buffer
//     const auto makeBatchBuffer = [] { return mappings_buffer(); };


//     //updates buffer with the database answer of a single query
//     const auto processQuery = [&] (mappings_buffer& buf,
//         sequence_query&& query, match_locations& allhits)
//     {
//         if(query.empty()) return;

//         if(opt.evaluate.determineGroundTruth) {
//             buf.truth.emplace_back(ground_truth_taxon(db, query.header));
//         }

//         // find the mapping for the query and add all results to coverage
//         buf.mappings.emplace_back();
//         for_each_eligible_mapping(db, opt.classify, query, allhits, [&](const auto& cand) {
//             buf.mappings.back().emplace_back(cand.tax);
//             buf.coverage.insert(allhits, cand);
//         });
//     };

//     //runs before a batch buffer is discarded
//     const auto finalizeBatch = [&] (mappings_buffer&& buf) {
//         if(opt.evaluate.determineGroundTruth) {
//             truth_.insert(truth_.end(), buf.truth.begin(), buf.truth.end());
//         }

//         mappings_.insert(mappings_.end(),
//                           std::make_move_iterator(buf.mappings.begin()),
//                           std::make_move_iterator(buf.mappings.end() ));

//         coverage_.merge(std::move(buf.coverage));
//     };

//     //runs if something needs to be appended to the output
//     const auto appendToOutput = [&] (const std::string&) {};

//     //run (parallel) database queries according to processing options
//     query_database(infiles, db, opt.process,
//                    makeBatchBuffer, processQuery, finalizeBatch,
//                    appendToOutput);

// //    std::cerr << "COVERAGEMIN: " << opt.classify.tgtCoverageMin << '\n';
//     coverage_filter(mappings_, coverage_, opt);

//     cout << '\n'
//          << "Hit Threshold: " << opt.classify.hitsMin << '\n'
//          << "Hit Cutoff: " << opt.classify.hitsDiffFraction << '\n'
//          << "Coverage Cutoff: " << opt.classify.tgtCoverageMin << '\n'
//          << "Coverage Max-Norm: " << int(opt.classify.tgtCoverageNorm) << endl;

//     evaluate_classification(opt, mappings_, truth_, results.statistics);
// }


// parameter search !
//
// they don't call it high *readability* computing
//
void map_queries_to_targets_rna_search_2pass(
    const vector<string>& infiles,
    const database& db, const query_options& opt)
{    
    struct eval_result {
        // parameters
        size_t hit_thresh = 0;
        double cutoff = 0;
        double cov_stat = 0;
        bool cov_maxnorm = false;

        // results
        size_t reads_total = 0;

        size_t matches_total = 0; 
        size_t reads_aligned = 0;
        size_t origin_mapped = 0;
        size_t correctly_rejected = 0;
    };  
    struct cmp {
        bool operator() (const eval_result& lhs, const eval_result& rhs) const noexcept{
            if (lhs.matches_total == rhs.matches_total)
                return lhs.origin_mapped > rhs.origin_mapped;
            else
                return lhs.matches_total < rhs.matches_total;
        }
    };
    std::multiset<eval_result, cmp> global_results;
    std::atomic<size_t> query_counter{0};
    matches_per_target_param pcoverage_;
    std::vector<std::pair<uint16_t, double>> settings;
    
    for (int rel = 50; rel <=100; rel+=2)
        for (int abs = 1; abs <=32; abs*=2)
            settings.emplace_back(abs, double(rel)/100);

    constexpr size_t N = 100;
    
    std::unordered_map<size_t, std::atomic<size_t>> total_matches, origin_mapped, aligned, correctly_rejected;
    std::unordered_map<size_t, std::atomic<size_t>> total_matches_maxnorm, origin_mapped_maxnorm, aligned_maxnorm, correctly_rejected_maxnorm;
    for (size_t i = 0; i < settings.size()*(N+1); ++i) {
        total_matches_maxnorm.emplace(i,0);
        origin_mapped_maxnorm.emplace(i,0);
        aligned_maxnorm.emplace(i,0);
        correctly_rejected_maxnorm.emplace(i,0);

        total_matches.emplace(i,0);
        origin_mapped.emplace(i,0);
        aligned.emplace(i,0);
        correctly_rejected.emplace(i,0);
    }

    std::vector<std::unordered_map<const taxon*, std::atomic<double>>> cov_cache;
    for (size_t i = 0; i < settings.size(); ++i) {
        cov_cache.emplace_back();
        for (const taxon& tax: db.target_taxa())
            cov_cache.back().emplace(&tax, -1.0);
    }

    const auto makeBatchBuffer = [] { return mappings_buffer(); };

    const auto processCoverage = [&] (mappings_buffer& buf,
        const sequence_query& query, const auto& allhits)
    {
        if(query.empty()) return;
        ++query_counter;

        buf.pcoverage.insert(allhits, make_classification_candidates(db, opt.classify, query, allhits));
    };

    //runs before a batch buffer is discarded
    const auto finalizeCoverage = [&] (mappings_buffer&& buf) {
        pcoverage_.merge(std::move(buf.pcoverage));
        std::cerr << query_counter << std::endl;
    };

    //runs if something needs to be appended to the output
    const auto appendToOutput = [&] (const std::string&) {
    };

    std::cerr << "STARTING FIRST PASS" << std::endl;
    //run (parallel) database queries according to processing options
    query_database(infiles, db, opt.pairing, opt.performance,
                   makeBatchBuffer, processCoverage, finalizeCoverage,
                   appendToOutput);
    
    query_counter = 0;

    const auto processQuery = [&] (mappings_buffer&,
        const sequence_query& query, const auto& allhits)
    {
        if(query.empty()) return;
        ++query_counter;
                
        classification cls { make_classification_candidates(db, opt.classify, query, allhits) } ;

        cls.groundTruth = ground_truth_taxon(db, query.header);

        for (size_t queue_ = 0; queue_ < settings.size(); ++queue_) {
            uint16_t abs = settings[queue_].first;
            double rel = settings[queue_].second;

            if (cls.candidates.empty()) {
                if (cls.groundTruth == nullptr) {
                    correctly_rejected_maxnorm[queue_*(N+1) + 0]++;
                    correctly_rejected[queue_*(N+1) + 0]++;
                }
                continue;
            }
            
            uint_least64_t max_hits = 0;
            for (const auto& cand: cls.candidates) max_hits = std::max(cand.hits, max_hits);

            double max_cov = 0;
            double min_cov = 1.0;
            for (const auto& cand: cls.candidates) {
                if (cand.hits < abs || double(cand.hits)/max_hits < rel)
                    continue;

                if (cov_cache[queue_][cand.tax] == -1.0)
                    cov_cache[queue_][cand.tax] = double(pcoverage_.num_hits(cand.tgt, abs, rel))/cand.tax->source().windows;
                max_cov = std::max(max_cov, cov_cache[queue_][cand.tax].load());
                min_cov = std::min(min_cov, cov_cache[queue_][cand.tax].load());
            }

            double max_cov_stat = 0;
            double min_cov_stat = 1.0;
            for (const auto& cand: cls.candidates) {
                if (cand.hits < abs || double(cand.hits)/max_hits < rel)
                    continue;
                
                double cov_stat = 1.0;
                if (max_cov > 0.0) {
                    cov_stat = cov_cache[queue_][cand.tax]/max_cov;
                }

                max_cov_stat = std::max(max_cov_stat, cov_stat);
                min_cov_stat = std::min(min_cov_stat, cov_stat);
                
                total_matches_maxnorm[cov_stat*N]++;
                total_matches[cov_cache[queue_][cand.tax]*N]++; 

                if (cls.groundTruth == cand.tax) {
                    origin_mapped_maxnorm[queue_*(N+1) + cov_stat*N]++;
                    origin_mapped[queue_*(N+1) + cov_cache[queue_][cand.tax]*N]++;
                }

            }

            aligned_maxnorm[queue_*(N+1) + max_cov_stat*N]++;
            aligned[queue_*(N+1) + max_cov*N]++;

            if (cls.groundTruth == nullptr) {
                correctly_rejected_maxnorm[queue_*(N+1) + std::ceil(min_cov_stat*N)]++;
                correctly_rejected[queue_*(N+1) + std::ceil(min_cov*N)]++;
            }
        }
    };
    
    //runs before a batch buffer is discarded
    const auto finalizeBatch = [&] (mappings_buffer&&) {
        std::cerr << query_counter << std::endl;
    };

    query_database(infiles, db, opt.pairing, opt.performance,
                   makeBatchBuffer, processQuery, finalizeBatch,
                   appendToOutput);

    // prefix_sums
    for (size_t queue_ = 0; queue_ < settings.size(); ++queue_) {
        for (size_t i = N-1; i < N; --i) {
            total_matches_maxnorm[queue_*(N+1) + i] += total_matches_maxnorm[queue_*(N+1) + i+1];
            origin_mapped_maxnorm[queue_*(N+1) + i] += origin_mapped_maxnorm[queue_*(N+1) + i+1];
            aligned_maxnorm[queue_*(N+1) + i] += aligned_maxnorm[queue_*(N+1) + i+1];
            correctly_rejected_maxnorm[queue_*(N+1) + N-i] += correctly_rejected_maxnorm[queue_*(N+1) + N-i-1];

            total_matches[queue_*(N+1) + i] += total_matches[queue_*(N+1) + i+1];
            origin_mapped[queue_*(N+1) + i] += origin_mapped[queue_*(N+1) + i+1];
            aligned[queue_*(N+1) + i] += aligned[queue_*(N+1) + i+1];
            correctly_rejected[queue_*(N+1) + N-i] += correctly_rejected[queue_*(N+1) + N-i-1];
        }
        for (size_t i = 0; i < N+1; ++i) {
            // filter here if neccessary

            uint16_t abs = settings[queue_].first;
            double rel = settings[queue_].second;
            global_results.insert({abs, rel, double(i)/N, true, query_counter, total_matches_maxnorm[queue_*(N+1) + i], aligned_maxnorm[queue_*(N+1) + i], origin_mapped_maxnorm[queue_*(N+1) + i], correctly_rejected_maxnorm[queue_*(N+1) + i]});
            global_results.insert({abs, rel, double(i)/N, false, query_counter, total_matches[queue_*(N+1) + i], aligned[queue_*(N+1) + i], origin_mapped[queue_*(N+1) + i], correctly_rejected[queue_*(N+1) + i]});
        }
    }

    std::cerr << global_results.size() << " results" << std::endl;
    std::cout << "hit_threshold,hit_cutoff,relative_coverage,reads,matches,aligned,TP,TN" << std::endl;
    for (const auto& i: global_results) {
        std::cout << i.hit_thresh << "," << i.cutoff << "," << i.cov_stat << "," << i.cov_maxnorm << ","
                << i.reads_total << "," << i.matches_total << "," << i.reads_aligned << "," << i.origin_mapped << "," << i.correctly_rejected << std::endl;
    }
    std::cout << std::endl;
}

// parameter search !
//
// they don't call it high *readability* computing
//
// from classify_rna.cpp
void map_queries_to_targets_rna_search_2pass(
    const vector<string>& infiles,
    const database& db, const query_options& opt)
{    
    struct eval_result {
        // parameters
        size_t hit_thresh = 0;
        double cutoff = 0;
        double cov_stat = 0;
        bool cov_maxnorm = false;

        // results
        size_t reads_total = 0;

        size_t matches_total = 0; 
        size_t reads_aligned = 0;
        size_t origin_mapped = 0;
        size_t correctly_rejected = 0;
    };  
    struct cmp {
        bool operator() (const eval_result& lhs, const eval_result& rhs) const noexcept{
            if (lhs.matches_total == rhs.matches_total)
                return lhs.origin_mapped > rhs.origin_mapped;
            else
                return lhs.matches_total < rhs.matches_total;
        }
    };
    std::multiset<eval_result, cmp> global_results;
    std::atomic<size_t> query_counter{0};
    matches_per_target_param pcoverage_;

    std::vector<std::pair<uint16_t, double>> settings;
    
    for (int rel = 0; rel <=100; rel+=1)
        for (int abs = 1; abs <=64; abs*=2)
            settings.emplace_back(abs, double(rel)/100);
    // settings.emplace_back(opt.classify.hitsMin, opt.classify.hitsCutoff);

    
    constexpr size_t N = 100;

    std::unordered_map<size_t, std::atomic<size_t>> total_matches, origin_mapped, aligned, correctly_rejected;
    std::unordered_map<size_t, std::atomic<size_t>> total_matches_maxnorm, origin_mapped_maxnorm, aligned_maxnorm, correctly_rejected_maxnorm;
    for (size_t i = 0; i < settings.size()*(N+1); ++i) {
        total_matches_maxnorm.emplace(i,0);
        origin_mapped_maxnorm.emplace(i,0);

        total_matches.emplace(i,0);
        origin_mapped.emplace(i,0);
        aligned.emplace(i,0);
        correctly_rejected.emplace(i,0);
    }

    for (size_t i = 0; i < settings.size(); ++i) { 
        correctly_rejected_maxnorm.emplace(i,0);
        aligned_maxnorm.emplace(i,0);
    }



    std::vector<std::unordered_map<const taxon*, std::atomic<double>>> cov_cache;
    for (size_t i = 0; i < settings.size(); ++i) {
        cov_cache.emplace_back();
        for (const taxon& tax: db.target_taxa())
            cov_cache.back().emplace(&tax, -1.0);
    }

    const auto makeBatchBuffer = [] { return mappings_buffer(); };

    const auto processCoverage = [&] (mappings_buffer& buf,
        const sequence_query& query, const auto& allhits)
    {
        if(query.empty()) return;
        ++query_counter;

        buf.pcoverage.insert(allhits, make_classification_candidates(db, opt.classify, query, allhits), opt.classify.covFill);
    };

    //runs before a batch buffer is discarded
    const auto finalizeCoverage = [&] (mappings_buffer&& buf) {
        pcoverage_.merge(std::move(buf.pcoverage));
        std::cerr << query_counter << std::endl;
    };

    //runs if something needs to be appended to the output
    const auto appendToOutput = [&] (const std::string&) {
    };

    std::cerr << "STARTING FIRST PASS" << std::endl;
    //run (parallel) database queries according to processing options
    query_database(infiles, db, opt.pairing, opt.performance,
                   makeBatchBuffer, processCoverage, finalizeCoverage,
                   appendToOutput);
    
    
    query_counter = 0;

    const auto processQuery = [&] (mappings_buffer&,
        const sequence_query& query, const auto& allhits)
    {
        if(query.empty()) return;
        ++query_counter;
                
        classification cls { make_classification_candidates(db, opt.classify, query, allhits) } ;
        cls.groundTruth = ground_truth_taxon(db, query.header);

        for (size_t queue_ = 0; queue_ < settings.size(); ++queue_) {
            
            uint16_t abs = settings[queue_].first;
            double rel = settings[queue_].second;

            auto cands = cls.candidates;
            auto filter_ops = opt.classify;
            filter_ops.hitsMin = abs;
            filter_ops.hitsCutoff = rel;
            hits_cutoff_filter(filter_ops, cands);

            size_t idx = (N+1) * queue_;
            
            if (cands.empty()) {
                if (cls.groundTruth == nullptr) {
                    correctly_rejected_maxnorm[queue_]++;
                    correctly_rejected[idx]++;
                }
                continue;
            }
            
            uint_least64_t max_hits = 0;
            for (const auto& cand: cands) max_hits = std::max(cand.hits, max_hits);

            double max_cov = 0;
            for (const auto& cand: cands) {
                if (cov_cache[queue_][cand.tax] == -1.0)
                    cov_cache[queue_][cand.tax] = double(pcoverage_.num_hits(cand.tgt, abs, rel))/cand.tax->source().windows;
                max_cov = std::max(max_cov, cov_cache[queue_][cand.tax].load());
            }

            for (const auto& cand: cands) {
                
                double cov_stat = cov_cache[queue_][cand.tax]/max_cov;

                total_matches_maxnorm[idx+N* cov_stat]++;
                total_matches[idx+N* cov_cache[queue_][cand.tax]]++; 

                if (cls.groundTruth == cand.tax) {
                    origin_mapped_maxnorm[idx+N* cov_stat]++;
                    origin_mapped[idx+N* cov_cache[queue_][cand.tax]]++;
                }
            }
            
            aligned_maxnorm[queue_]++;
            aligned[idx+N* max_cov]++;

            if (cls.groundTruth == nullptr) {
                size_t reject_idx = idx+N*max_cov +1;
                if (reject_idx <= idx+N)
                    correctly_rejected[reject_idx]++;
            }
        }
    };

    //runs before a batch buffer is discarded
    const auto finalizeBatch = [&] (mappings_buffer&&) {
        std::cerr << query_counter << std::endl;
    };

    query_database(infiles, db, opt.pairing, opt.performance,
                   makeBatchBuffer, processQuery, finalizeBatch,
                   appendToOutput);

    // prefix_sums

    do_in_parallel(opt.performance.numThreads, settings.size(), [&](size_t q) {
        size_t idx = (N+1) * q;
        for (size_t i = N-1; i < N; --i) {
            total_matches_maxnorm.at(idx+ i) += total_matches_maxnorm.at(idx+ i+1);
            origin_mapped_maxnorm[idx+ i] += origin_mapped_maxnorm[idx+ i+1];

            total_matches[idx+ i] += total_matches[idx+ i+1];
            origin_mapped[idx+ i] += origin_mapped[idx+ i+1];
            aligned[idx+ i] += aligned[idx+ i+1];
            correctly_rejected[idx+ N-i] += correctly_rejected[idx+ N-i-1];
        }
    });

    for (size_t q = 0; q < settings.size(); ++q) {
        size_t idx = (N+1) * q;
        for (size_t i = 0; i <= N; ++i) {
            // filter here if neccessary

            uint16_t abs = settings[q].first;
            double rel = settings[q].second;
            global_results.insert({abs, rel, double(i)/N, true, query_counter, total_matches_maxnorm[idx+ i], aligned_maxnorm[q], origin_mapped_maxnorm[idx+ i], correctly_rejected_maxnorm[q]});
            global_results.insert({abs, rel, double(i)/N, false, query_counter, total_matches[idx+ i], aligned[idx+ i], origin_mapped[idx+ i], correctly_rejected[idx+ i]});
        }
    }

    std::cerr << global_results.size() << " results" << std::endl;
    std::cout << "hit_threshold,hit_cutoff,relative_coverage,reads,matches,aligned,TP,TN" << std::endl;
    for (const auto& i: global_results) {
        std::cout << i.hit_thresh << "," << i.cutoff << "," << i.cov_stat << "," << i.cov_maxnorm << ","
                << i.reads_total << "," << i.matches_total << "," << i.reads_aligned << "," << i.origin_mapped << "," << i.correctly_rejected << std::endl;
    }
    std::cout << std::endl;
}
