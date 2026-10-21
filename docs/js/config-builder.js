document.addEventListener("DOMContentLoaded", () => {

    const downloadButton =
        document.getElementById("download-config");

    downloadButton.addEventListener(
        "click",
        generateConfig
    );

});

function generateConfig() {

    const project =
        document.getElementById("projectName").value;

    const forward =
        document.getElementById("forwardReads").value;

    const reverse =
        document.getElementById("reverseReads").value;

    const threads =
        document.getElementById("threads").value;

    const quality =
        document.getElementById("quality").value;

    const config = `
# MJOLNIR3 Configuration

PROJECT=${project}

FORWARD_READS=${forward}

REVERSE_READS=${reverse}

THREADS=${threads}

QUALITY=${quality}
`;

    const blob =
        new Blob([config], {
            type: "text/plain"
        });

    const link =
        document.createElement("a");

    link.href =
        URL.createObjectURL(blob);

    link.download =
        "mjolnir3_config.txt";

    link.click();

    URL.revokeObjectURL(link.href);
}