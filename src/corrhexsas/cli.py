"""Command-line entry point; all computation is local."""
import argparse
import json
import sys
from . import __version__
from .examples import create_example
from .workflow import run_fit, check_convergence


def main(argv=None):
    parser = argparse.ArgumentParser(description="Fit hexagonal SAS using correlated displacements and finite grains.")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    ex = sub.add_parser("example", help="Create synthetic data and a runnable fitting configuration")
    ex.add_argument("--output", required=True)
    fit = sub.add_parser("fit", help="Fit one or more q I dI [sigma_q] files")
    fit.add_argument("config")
    fit.add_argument("--output", required=True)
    fit.add_argument("--plot", action="store_true")
    check = sub.add_parser("check", help="Compare a saved fit with refined numerical settings, without refitting")
    check.add_argument("config")
    check.add_argument("--fit", required=True, dest="fit_result")
    check.add_argument("--output", required=True)
    check.add_argument("--factor", type=float, default=1.5)
    args = parser.parse_args(argv)
    try:
        if args.command == "example":
            print(create_example(args.output))
            return 0
        if args.command == "fit":
            result = run_fit(args.config, args.output, plot=args.plot)
            print(json.dumps({k: result[k] for k in ("success", "objective_name", "objective", "chi2", "parameters")}, indent=2))
            return 0 if result["success"] else 2
        result = check_convergence(args.config, args.fit_result, args.output, args.factor)
        print(json.dumps(result["metrics"], indent=2))
        return 0
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        print(f"corrhexsas: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
