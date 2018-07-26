import re


class Snippet:
    def __init__(self, tag, start, end):
        self.tag = tag
        self.start = start
        self.end = end

    def matches_start(self, string):
        return self.start is not None and string.startswith(self.start)

    def matches_end(self, string):
        return self.end is not None and string.startswith(self.end)

    def get_tag(self):
        return "/**! [Snippet%s] */" % self.tag


def handle_line(line, snippet, target_file):
    start_pattern_matches = snippet.matches_start(line)
    end_pattern_matches = snippet.matches_end(line)

    if start_pattern_matches or end_pattern_matches:
        target_file.write("%s\n" % snippet.get_tag())

    return end_pattern_matches


none_snippet = Snippet("", None, None)

snippet_list = [
    Snippet("Version", "SCIP version", "SCIP> help"),
    Snippet("Help", "SCIP> help", "SCIP> read check/instances/MIP/stein27.fzn"),
    Snippet("Opt1", "SCIP> read check/instances/MIP/stein27.fzn", "SCIP> write solution stein27.sol"),
    Snippet("WriteSolutions", "SCIP> write solution stein27.sol", "SCIP> set save settingsfile.set"),
    Snippet("SaveSettingsOverview", "SCIP> set save settingsfile.set", "SCIP> display heuristics"),
    Snippet("DisplayStatistics", "SCIP> display heuristics", "SCIP> set"),
    Snippet("SetSettings", "SCIP> set", "SCIP> set default"),
    Snippet("Opt2", "SCIP> set default", "SCIP> set save settingsfile.set"),
    Snippet("SaveSettingsFull", "SCIP> set save settingsfile.set", "SCIP> set diffsave settingsfile.set"),
    Snippet("SaveSettingsDiff", "SCIP> set diffsave settingsfile.set", "SCIP> set load settingsfile.set"),
    Snippet("LoadSettings", "SCIP> set load settingsfile.set", "SCIP> quit")
]


def main():
    with open("inc/shelltutorial/shelltutorialannotated.tmp", "w") as annotated_file:
        with open("inc/shelltutorial/shelltutorialraw.tmp", "r") as raw_file:
            snippet_iterator = iter(snippet_list)
            current_snippet = next(snippet_iterator)
            for line in raw_file:

                while handle_line(line, current_snippet, annotated_file):
                    try:
                        current_snippet = next(snippet_iterator)
                    except StopIteration:
                        current_snippet = none_snippet

                annotated_file.write(line)


if __name__ == "__main__":
    main()
