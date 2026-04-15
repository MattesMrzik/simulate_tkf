use anyhow::Context;
use clap::Parser;
use memory_stats::memory_stats;
use phylo::Result;
use phylo::alignment::{Alignment, AlignmentSimulation, AncestralAlignment, MASA};
use phylo::io::{read_newick_from_file, write_newick_to_file};
use phylo::random::DefaultGenerator;
use phylo::substitution_models::{JC69, SubstModel};
use phylo::tkf_model::TKF92IndelModel;
use phylo::tkf_model::sim_tkf_msa::TKFMSASimulator;
use phylo::tree::Tree;
use std::collections::HashSet;
use std::fs;
use std::io::Write;

mod args;
use crate::args::Args;

pub(crate) fn set_missing_tree_node_ids(tree: &Tree) -> Result<Tree> {
    let mut tree_with_all_ids = tree.clone();
    let mut seen_user_set_ids = HashSet::new();
    let mut count = 0;
    for node_idx in tree.postorder() {
        let id = tree.node_id(node_idx);
        if id.is_empty() {
            let mut new_id = format!("I{count}");
            while !seen_user_set_ids.insert(new_id.clone()) {
                count += 1;
                new_id = format!("I{count}");
            }
            tree_with_all_ids.node_mut(node_idx).id = new_id.clone();
        }
    }
    Ok(tree_with_all_ids)
}

fn write_info_file(
    path: &std::path::Path,
    msa_len: usize,
    time_ms: u128,
    memory_mb: f64,
    seed: u64,
    args: &Args,
    attempts: Option<u32>,
) {
    let mut file = fs::File::create(path).expect("Unable to create info file");
    writeln!(file, "msa_len: {}", msa_len).unwrap();
    writeln!(file, "time_ms: {}", time_ms).unwrap();
    writeln!(file, "memory_mb: {:.4}", memory_mb).unwrap();
    writeln!(file, "seed: {}", seed).unwrap();
    writeln!(file, "lambda: {}", args.lambda).unwrap();
    writeln!(file, "mu: {}", args.mu).unwrap();
    writeln!(file, "r: {}", args.r).unwrap();
    if let Some(root_len) = args.root_length {
        writeln!(file, "root_length: {}", root_len).unwrap();
    }
    writeln!(
        file,
        "tree_file: {}",
        args.tree_file.to_str().unwrap_or("invalid_path")
    )
    .unwrap();
    if let Some(attempts) = attempts {
        writeln!(file, "attempts: {}", attempts).unwrap();
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    let newick = fs::read_to_string(&args.tree_file).expect("Unable to read tree file");
    let tree = phylo::tree::tree_parser::from_newick(&newick)
        .expect("Unable to parse Newick tree")
        .pop()
        .expect("Tree file was empty");
    let tree = set_missing_tree_node_ids(&tree).unwrap();

    let subst_model = SubstModel::<JC69>::new(&[], &[]);

    let tkf_model = TKF92IndelModel::new(args.lambda, args.mu, args.r);

    let rng = match args.seed {
        Some(seed) => DefaultGenerator::new(seed),
        None => DefaultGenerator::default(),
    };
    let seed_used = rng.seed();

    let mut simulator = TKFMSASimulator::new(
        tkf_model,
        subst_model,
        tree.clone(),
        rng,
        args.max_insertion_length,
    );
    if let Some(root_len) = args.root_length {
        simulator.root_length(Some(root_len));
    }

    fs::create_dir_all(&args.output_dir).expect("Unable to create output directory");

    // Simulate until we get a valid MSA (length > 0 and no all-gap leaf sequences)
    let mut initial_mem;
    let mut final_mem;
    let mut start;
    let mut duration;
    let mut msa;
    let mut attempt = 0;

    loop {
        attempt += 1;
        initial_mem = memory_stats().map(|ms| ms.physical_mem).unwrap_or(0);
        start = std::time::Instant::now();
        msa = simulator.simulate_ancestral_alignment::<MASA>();
        duration = start.elapsed();
        final_mem = memory_stats().map(|ms| ms.physical_mem).unwrap_or(0);

        let mem_diff = final_mem as i64 - initial_mem as i64;
        let mem_mb = mem_diff as f64 / 1024.0 / 1024.0;

        // Write info for this attempt
        write_info_file(
            &args.output_dir.join(format!("info_{}.txt", attempt)),
            msa.len(),
            duration.as_millis(),
            mem_mb,
            seed_used,
            &args,
            None,
        );

        if msa.len() == 0 {
            continue;
        }

        let mut all_leaves_have_chars = true;
        for leaf_record in msa.seqs() {
            if leaf_record.seq().is_empty() {
                all_leaves_have_chars = false;
                break;
            }
        }

        if all_leaves_have_chars {
            break;
        }
    }

    // Output basic info about the simulated MSA
    println!("Simulation took: {:?}", duration);
    println!("MSA length: {}", msa.len());

    let mem_diff = final_mem as i64 - initial_mem as i64;
    let mem_mb = mem_diff as f64 / 1024.0 / 1024.0;
    println!("Memory usage of simulation: {:.2} MB", mem_mb);

    // Write MASA to masa.fasta
    let mut masa_file =
        fs::File::create(args.output_dir.join("masa.fasta")).expect("Unable to create masa.fasta");
    write!(masa_file, "{}", msa).expect("Unable to write to masa.fasta");

    // Convert MASA to MSA (leaf nodes only) and write to msa.fasta
    let leaf_msa: phylo::alignment::MSA = msa.clone().into_alignment(&tree);
    let mut msa_file =
        fs::File::create(args.output_dir.join("msa.fasta")).expect("Unable to create msa.fasta");
    write!(msa_file, "{}", leaf_msa).expect("Unable to write to msa.fasta");

    // Final successful info.txt
    write_info_file(
        &args.output_dir.join("info.txt"),
        msa.len(),
        duration.as_millis(),
        mem_mb,
        seed_used,
        &args,
        Some(attempt),
    );

    for node in tree.preorder() {
        println!("Node ID: {}", tree.node_id(node),);
    }

    write_newick_to_file(
        std::slice::from_ref(&tree),
        args.output_dir.join("tree.nwk"),
    )
    .context("Failed to write optimized tree to file")?;
    let tree_from_just_written_file = read_newick_from_file(args.output_dir.join("tree.nwk"))
        .expect("Unable to read back the just written tree file")
        .pop()
        .expect("The just written tree file was empty");
    println!(
        "Tree read back from file matches original: {}",
        tree_from_just_written_file
    );
    println!("newick = {}", tree.to_newick());

    println!("Results written to: {:?}", args.output_dir);

    if args.ali_view {
        let msa_path = args.output_dir.join("masa.fasta");
        println!("Opening results in AliView...");
        std::process::Command::new("open")
            .args(["-g", "-a", "AliView"]) // -a is for specifying the application, -g is for opening in the background
            .arg(msa_path)
            .status()
            .expect("Failed to open AliView");
    }
    Ok(())
}
