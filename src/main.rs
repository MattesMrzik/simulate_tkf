#![allow(dead_code)]

use anyhow::{Context, Result, anyhow};
use clap::Parser;
use memory_stats::memory_stats;
use std::collections::HashSet;
use std::fs;
use std::io::Write;

use phylo::alignment::{Alignment, AlignmentSimulation, AncestralAlignment, MASA};
use phylo::io::write_newick_to_file;
use phylo::random::DefaultGenerator;
use phylo::substitution_models::{JC69, SubstModel};
use phylo::tkf_model::TKF92IndelModel;
use phylo::tkf_model::simulate_msa::{RootLength, TKFMSASimulator};
use phylo::tree::Tree;

use crate::args::Args;

mod args;

fn has_any_leaf(msa: &MASA) -> bool {
    for leaf_record in msa.seqs() {
        if !leaf_record.seq().is_empty() {
            return true;
        }
    }
    false
}

fn all_leaves_have_chars(msa: &MASA) -> bool {
    for leaf_record in msa.seqs() {
        if leaf_record.seq().is_empty() {
            return false;
        }
    }
    true
}

fn set_missing_tree_node_ids(tree: &Tree) -> Tree {
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
    tree_with_all_ids
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
    if let Some(root_len) = &args.root_length {
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
    let tree = set_missing_tree_node_ids(&tree);

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
    if let Some(root_len) = args.root_length.as_ref() {
        let s = root_len.to_string();
        let root_len = match s.to_lowercase().as_str() {
            "sampled" => RootLength::Sampled,
            "expected" => RootLength::Expected,
            _ => {
                let n = s.parse::<usize>().map_err(|_| {
                    anyhow!(
                        "invalid root length '{}'; expected 'sampled', 'expected', or a positive integer",
                        s
                    )
                })?;
                RootLength::Defined(n)
            }
        };
        simulator.root_length(root_len);
    }

    fs::create_dir_all(&args.output_dir).expect("Unable to create output directory");

    let mut masa;
    let mut duration;
    let mut attempt = 0;

    loop {
        let initial_mem = memory_stats().map(|ms| ms.physical_mem).unwrap_or(0);
        let start = std::time::Instant::now();
        masa = simulator.simulate_ancestral_alignment::<MASA>();
        duration = start.elapsed();
        let final_mem = memory_stats().map(|ms| ms.physical_mem).unwrap_or(0);

        if masa.len() == 0 {
            continue;
        }
        attempt += 1;

        let mem_diff = final_mem as i64 - initial_mem as i64;
        let mem_mb = mem_diff as f64 / 1024.0 / 1024.0;

        write_info_file(
            &args.output_dir.join(format!("info_{}.txt", attempt)),
            masa.len(),
            duration.as_millis(),
            mem_mb,
            seed_used,
            &args,
            None,
        );

        if all_leaves_have_chars(&masa) || !args.retry_if_empty_leaf {
            break;
        }
    }

    println!("Simulation took: {:?}", duration);
    println!("MSA length: {}", masa.len());
    masa.remove_extinct_columns();

    let masa_len = masa.len();
    let mut masa_file =
        fs::File::create(args.output_dir.join("masa.fasta")).expect("Unable to create masa.fasta");
    write!(masa_file, "{}", masa).expect("Unable to write to masa.fasta");

    let leaf_msa: phylo::alignment::MSA = masa.clone().into_alignment(&tree);
    let leaf_msa_len = leaf_msa.len();
    assert!(
        leaf_msa_len == masa_len,
        "Leaf MSA length should match ancestral MSA length, since i called remove_extinct_columns() on the ancestral MSA"
    );
    let mut msa_file =
        fs::File::create(args.output_dir.join("msa.fasta")).expect("Unable to create msa.fasta");
    write!(msa_file, "{}", leaf_msa).expect("Unable to write to msa.fasta");

    // copy the latest info file to info.txt for easy access
    fs::copy(
        args.output_dir.join(format!("info_{}.txt", attempt)),
        args.output_dir.join("info.txt"),
    )
    .expect("Unable to copy info file");

    write_newick_to_file(
        std::slice::from_ref(&tree),
        args.output_dir.join("tree.nwk"),
    )
    .context("Failed to write tree to file")?;

    // TODO: also write the other newick with and without brackets around the root, see root_tree

    println!("Results written to: {:?}", args.output_dir);

    if args.ali_view {
        let msa_path = args.output_dir.join("masa.fasta");
        println!("Opening results in AliView...");
        std::process::Command::new("open")
            .args(["-g", "-a", "AliView"])
            .arg(msa_path)
            .status()
            .expect("Failed to open AliView");
    }
    Ok(())
}
