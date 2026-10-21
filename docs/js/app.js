document.addEventListener("DOMContentLoaded", () => {

    const helpButton =
        document.getElementById("toggle-help");

    const helpSection =
        document.getElementById("help-section");

    helpButton.addEventListener("click", () => {

        helpSection.classList.toggle("active");

        if (helpSection.classList.contains("active")) {

            helpButton.textContent = "Close Help";

        } else {

            helpButton.textContent = "Open Help";
        }
    });

});