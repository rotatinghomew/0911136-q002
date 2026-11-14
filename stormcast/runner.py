"""运行领域摘要、极端降雨频率分析或启动本地供数服务。"""

import argparse
import json
from pathlib import Path
from http.server import ThreadingHTTPServer

from .analysis import summarize
from .datasource import SourceClient
from .frequency import frequency_summary
from .reports import package_report, save_report
from .service import handler_for, load_repository, serve


def build_summary(tables, mode, bootstrap_replicates):
    """mode: 'legacy' 仅现行经验摘要；'frequency' 仅频率分析；'all' 两者合并。"""
    legacy = summarize(tables)
    if mode == "legacy":
        return legacy
    frequency = frequency_summary(tables, replicates=bootstrap_replicates)
    if mode == "frequency":
        return frequency
    return {"empirical_summary": legacy, "frequency_analysis": frequency}


def main():
    parser = argparse.ArgumentParser(description="领域资料与极端降雨频率分析")
    parser.add_argument("--dataset", default="fixtures/dataset.json")
    parser.add_argument("--source-url")
    parser.add_argument("--output", default="/tmp/domain-summary.json")
    parser.add_argument("--serve", action="store_true")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--mode", choices=("all", "legacy", "frequency"),
                        default="all",
                        help="all=现行经验摘要加频率分析（默认）；legacy=仅现行摘要；frequency=仅频率分析")
    parser.add_argument("--bootstrap-replicates", type=int, default=1000,
                        help="尾部重现期自助区间的重采样次数")
    arguments = parser.parse_args()
    schema_path = Path(__file__).with_name("schema.json")
    if arguments.serve:
        server = ThreadingHTTPServer(("127.0.0.1", arguments.port), handler_for(load_repository(arguments.dataset, schema_path)))
        try:
            print(json.dumps({"listening": server.server_address}, ensure_ascii=False), flush=True)
            server.serve_forever()
        finally:
            server.server_close()
        return
    if arguments.source_url:
        snapshot = SourceClient(arguments.source_url).read()
    else:
        with serve(load_repository(arguments.dataset, schema_path)) as base_url:
            snapshot = SourceClient(base_url).read()
    summary = build_summary(snapshot.tables, arguments.mode,
                            arguments.bootstrap_replicates)
    report = package_report(snapshot, summary)
    save_report(arguments.output, report)
    print(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
