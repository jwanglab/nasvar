use std::time::{SystemTime, UNIX_EPOCH};

/// Get current UTC timestamp in ISO 8601 format (e.g. "2025-02-05T14:30:00Z")
pub fn utc_now_iso8601() -> String {
    let (year, month, day, hours, minutes, seconds) = utc_now_parts();
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        year, month, day, hours, minutes, seconds
    )
}

fn utc_now_parts() -> (i64, u32, u32, u64, u64, u64) {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    let secs_per_day = 86400u64;
    let days = secs / secs_per_day;
    let secs_today = secs % secs_per_day;

    let hours = secs_today / 3600;
    let minutes = (secs_today % 3600) / 60;
    let seconds = secs_today % 60;

    let (year, month, day) = days_to_ymd(days as i64);
    (year, month, day, hours, minutes, seconds)
}

fn days_to_ymd(days: i64) -> (i64, u32, u32) {
    let mut remaining = days;
    let mut year = 1970i64;

    loop {
        let days_in_year = if is_leap_year(year) { 366 } else { 365 };
        if remaining < days_in_year {
            break;
        }
        remaining -= days_in_year;
        year += 1;
    }

    let leap = is_leap_year(year);
    let days_in_month: [i64; 12] = if leap {
        [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    } else {
        [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    };

    let mut month = 1u32;
    for dim in days_in_month.iter() {
        if remaining < *dim {
            break;
        }
        remaining -= *dim;
        month += 1;
    }

    let day = remaining as u32 + 1;
    (year, month, day)
}

fn is_leap_year(year: i64) -> bool {
    (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iso8601_format() {
        let ts = utc_now_iso8601();
        assert!(ts.contains('T'));
        assert!(ts.ends_with('Z'));
        assert_eq!(ts.len(), 20);
    }

}
