

from pathlib import Path
from ipyrad2.utils.parallel import run_with_pool, run_pipeline
from loguru import logger



def _test_func(a, b):
    # time.sleep(2)
    raise ValueError("x")
    return a + b


if __name__ == "__main__":

    jobs = {
        "A": (run_pipeline, {"cmds": [["ls", "-l"], ["wc"]]}),
        "B": (run_pipeline, {"cmds": [["ls", "-l"]], "outfile": Path("/tmp/test.txt")}),
    }
    res = run_with_pool(jobs, "INFO", 4)
    print(res)



    # jobs = {
    #     "A": {"a": 3, "b": 4},
    #     "B": {"a": 3, "b": 8},
    # }
    # try:
    #     res = run_with_pool(_test_func, jobs, 4)
    #     print(res)
    # except KeyboardInterrupt:
    #     print('stopped gracefully')
    # except ValueError as exc:
    #     logger.error(exc)
    # except Exception as exc:
    #     logger.exception(exc)

    # print(run_pipeline([["ls", "-l"], ["wc"]]))

    # print(run_pipeline([["echo", "hello world"]]))
