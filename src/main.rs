use clap::Parser;
use rust_htslib::bam;

pub mod align_stats;

#[derive(Parser, Debug)]
#[clap(author, version, about="Quickly compute long read stats for alignment files from stdin", long_about = None)]
struct Cli {
    /// out prefix string to use for naming output files
    #[arg(short, long, value_parser, default_value_t = String::from("out"))]
    out_prefix: String,
}

fn main() {
    let args = Cli::parse();
    if atty::is(atty::Stream::Stdin) {
        use clap::CommandFactory;
        Cli::command().print_help().unwrap();
        eprintln!("\nCould not open bam from stdin. Pipe samtools view output from stdin!\n");
        std::process::exit(1);
    }

    run_lra_stats(&args);
}

fn run_lra_stats(args: &Cli) {
    let mut bam = bam::Reader::from_stdin().expect("could not open bam/cram from stdin");
    let full_stats = align_stats::extract(&mut bam);
    let agg_stats = full_stats.aggregate();
    agg_stats.write_as_json(args.out_prefix.clone() + ".lra_stats_aggregate.json");
    full_stats.write_as_arrow(args.out_prefix.clone() + ".lra_stats_nanoplot.feather");
}
