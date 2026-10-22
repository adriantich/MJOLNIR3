// document.addEventListener("DOMContentLoaded", () => {

//     const buttons =
//         document.querySelectorAll(".help-btn");

//     buttons.forEach(button => {

//         button.addEventListener("click", () => {

//             const targetId =
//                 button.dataset.target;

//             const helpPanel =
//                 document.getElementById(targetId);

//             helpPanel.classList.toggle("active");

//             button.textContent =
//                 helpPanel.classList.contains("active")
//                 ? "Close Help"
//                 : "Open Help";

//         });

//     });

// });

document
    .querySelectorAll(".help-btn")
    .forEach(button => {

        button.addEventListener("click", () => {

            const formGroup =
                button.closest(".form-group");

            const helpSection =
                formGroup.querySelector(
                    ".help-section"
                );

            helpSection.classList.toggle(
                "active"
            );

        });

    });

// document
//     .querySelectorAll(".help-btn")
//     .forEach(button => {

//         button.addEventListener("click", () => {

//             const stepHeader =
//                 button.closest(".step-header");

//             const helpSectionHeader =
//                 stepHeader.querySelector(
//                     ".help-section"
//                 );

//             helpSectionHeader.classList.toggle(
//                 "active"
//             );

//         });

//     });