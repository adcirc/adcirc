#!/usr/bin/env python3

# This file implements a series of "fuzz" tests which
# validate that the fdate library produces identical results
# to the python datetime library by generating random date information
# and performing some work on it (i.e., parsing or date arithmetic)

from datetime import datetime, timezone, timedelta
import random

random.seed()


def random_date() -> datetime:
    random_month = random.randint(1, 12)
    if random_month == 2:
        random_day = random.randint(1, 28)
    elif random_month in [4, 6, 9, 11]:
        random_day = random.randint(1, 30)
    else:
        random_day = random.randint(1, 31)

    return datetime(
        year=random.randint(1900, 2100),
        month=random_month,
        day=random_day,
        hour=random.randint(0, 23),
        minute=random.randint(0, 59),
        second=random.randint(0, 59),
        tzinfo=timezone.utc,
    )


def random_timedelta() -> timedelta:
    day_range = 1000
    hour_range = 1000
    minute_range = 1000
    second_range = 1000
    return timedelta(
        days=random.randint(-day_range, day_range),
        hours=random.randint(-hour_range, hour_range),
        minutes=random.randint(-minute_range, minute_range),
        seconds=random.randint(-second_range, second_range),
    )


def run_arithmetic_tests(args):
    import subprocess

    # If the user has specified a timeout, then the number of tests is effectively infinite
    test_count = args.test_count if args.timeout is None else int(1e6)

    # To track the timeout, we get the current time and add a timeout to it
    if args.timeout is not None:
        timeout_time = datetime.now(tz=timezone.utc) + timedelta(seconds=args.timeout)
    else:
        timeout_time = datetime(2050, 1, 1, tzinfo=timezone.utc)

    for i in range(test_count):

        if args.timeout is not None:
            # Check if the current time has exceeded the timeout
            current_time = datetime.now(tz=timezone.utc)
            if current_time > timeout_time:
                print(f"[INFO] Test timeout ({args.timeout} seconds) reached after {i} tests.")
                break

        random_start_time = random_date()
        random_delta = random_timedelta()
        python_result = random_start_time + random_delta
        python_result_str = python_result.strftime("%Y-%m-%dT%H:%M:%S")

        date_str = random_start_time.strftime("%Y-%m-%dT%H:%M:%S")
        delta_str_days = random_delta.days
        delta_str_hours = random_delta.seconds // 3600
        delta_str_minutes = (random_delta.seconds - delta_str_hours * 3600) // 60
        delta_str_seconds = (
            random_delta.seconds - delta_str_hours * 3600 - delta_str_minutes * 60
        )

        try:
            call_args = [
                args.exe,
                "arithmetic",
                date_str,
                str(delta_str_days),
                str(delta_str_hours),
                str(delta_str_minutes),
                str(delta_str_seconds),
            ]
            result = subprocess.run(
                call_args,
                capture_output=True,
                text=True,
                check=True,
            )
            output_time_str = result.stdout.strip()

            if output_time_str != python_result_str:
                print(
                    f"[{i+1}/{args.test_count}] Test with start time: {date_str} and delta: {delta_str_days} days, {delta_str_hours} hours, {delta_str_minutes} minutes, {delta_str_seconds} seconds"
                )
                print(
                    f"    Mismatch: Expected [python]'{python_result_str}', got [c++]'{output_time_str}'"
                )
                raise RuntimeError(
                    f"Mismatch in output for start time: {date_str} and delta: {delta_str_days} days, {delta_str_hours} hours, {delta_str_minutes} minutes, {delta_str_seconds} seconds"
                )

        except subprocess.CalledProcessError as e:
            diag_str = " ".join(call_args)
            raise RuntimeError(
                f"Error running the executable: {args.exe} with {diag_str}"
            )


def run_parse_tests(args):
    import subprocess
    import random

    format_strings = [
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M:%SZ",
        "%Y-%m-%dT%H:%M:%S",
        "%Y/%m/%d %H:%M:%S",
        "%Y.%m.%d %H:%M:%S",
        "%Y%m%d%H%M%S",
        "%Y/%m/%d %H:%M",
        "%Y-%m-%d",
        "%Y/%m/%d",
        "%Y.%m.%d",
        "%Y%m%d",
    ]

    # If the user has specified a timeout, then the number of tests is effectively infinite
    test_count = args.test_count if args.timeout is None else int(1e6)

    # To track the timeout, we get the current time and add a timeout to it
    if args.timeout is not None:
        timeout_time = datetime.now(tz=timezone.utc) + timedelta(seconds=args.timeout)
    else:
        timeout_time = datetime(2050, 1, 1, tzinfo=timezone.utc)

    for i in range(test_count):

        if args.timeout is not None:
            # Check if the current time has exceeded the timeout
            current_time = datetime.now(tz=timezone.utc)
            if current_time > timeout_time:
                print(f"[INFO] Test timeout ({args.timeout} seconds) reached after {i} tests.")
                break

        random_datetime = random_date()
        format_str = random.choice(format_strings)
        datetime_str = random_datetime.strftime(format_str)
        parsed_datetime = datetime.strptime(datetime_str, format_str).replace(
            tzinfo=timezone.utc
        )
        timestamp = int(parsed_datetime.timestamp()) * 1000

        try:
            result = subprocess.run(
                [args.exe, "parse", datetime_str],
                capture_output=True,
                text=True,
                check=True,
            )
            output_timestamp = int(result.stdout.strip())
            if output_timestamp != timestamp:
                print(
                    f"[{i+1}/{args.test_count}] Test with datetime: {datetime_str} and format: {format_str}"
                )
                print(
                    f"    Mismatch: Expected [python]'{timestamp}', got [c++]'{output_timestamp}'"
                )
                raise RuntimeError(
                    f"Mismatch in output for datetime: {datetime_str} and format: {format_str}"
                )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Error running the executable: {args.exe} with datetime: {datetime_str} and format: {format_str}"
            ) from e

    print(f"All {args.test_count} tests passed successfully!")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Fuzz datetime interface")
    parser.add_argument(
        "--test-count",
        type=int,
        default=100,
        help="Number of tests to run (default: 100)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=None,
        help="Timeout for each test in seconds (default: None, no timeout)",
    )
    parser.add_argument(
        "--exe", type=str, required=True, help="Path to the C++ executable"
    )
    parser.add_argument(
        "--test",
        required=True,
        help="Type of test to run",
        choices=["parse", "arithmetic"],
    )
    args = parser.parse_args()

    if args.test == "parse":
        run_parse_tests(args)
    elif args.test == "arithmetic":
        run_arithmetic_tests(args)
    else:
        raise ValueError(
            f"Unknown test type: {args.test}. Supported types are 'parse' and 'arithmetic'."
        )


if __name__ == "__main__":
    main()
