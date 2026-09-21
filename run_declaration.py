"""Run the v3 declared-set robust model, with optional confirmed full-cycle feedback."""
import argparse
import json
from pathlib import Path

from lcdl_algorithm import Case, SolverOptions, run_two_stage, export_result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,required=True)
    parser.add_argument('--feedback',type=Path,help='UTF-8 JSON containing a K x T MW array')
    parser.add_argument('--output',type=Path,default=Path('results/declaration/result.json'))
    parser.add_argument('--time-limit',type=float,default=60)
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--skip-cost',action='store_true',help='Verify robust feasibility only; Eq. (31) is optional')
    args=parser.parse_args()
    case=Case.from_json(args.case)
    if not case.declaration_model: parser.error('Use schema_version=3; migrate rho/s explicitly')
    feedback=json.loads(args.feedback.read_text(encoding='utf-8')) if args.feedback else None
    result=run_two_stage(case,SolverOptions(time_limit=args.time_limit,threads=args.threads,
                                           optimize_robust_cost=not args.skip_cost),feedback_p_mw=feedback)
    output=export_result(result,args.output,include_matrices=True)
    print(result['status'],result['reason'])
    print(output)
    return 0 if result['declaration_feasible'] and (feedback is None or result['executable']) else 2


if __name__=='__main__':
    raise SystemExit(main())
