#include "raylib.h"
#include <stdlib.h>         
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

#define MAX_FILEPATH_SIZE 2048

#define SCREEN_FACTOR 100
#define S_WIDTH 16*SCREEN_FACTOR
#define S_HEIGHT 9*SCREEN_FACTOR
#define NUM_ROWS 12
#define ROW S_HEIGHT / NUM_ROWS

#define LARGE_TEXT 40
#define MEDIUM_TEXT 20
#define SMALL_TEXT 10


int main(void)
{
   
    InitWindow(S_WIDTH , S_HEIGHT, "raylib [core] example - drop files");

    bool file_is_dropped = false;
    bool multi_files_dropped = false;
    bool run_button = false;
    char file_path[MAX_FILEPATH_SIZE]; 

    SetTargetFPS(60);               
    
    while (!WindowShouldClose())    
    {
        if (IsFileDropped()) {
            FilePathList dropped_file = LoadDroppedFiles();
            if (dropped_file.count > 1) {
                multi_files_dropped = true;
            }
            else {
                TextCopy(file_path, dropped_file.paths[0]);
                file_is_dropped = true;
                multi_files_dropped = false;
            }
            UnloadDroppedFiles(dropped_file);    
        }

        BeginDrawing();
            ClearBackground(RAYWHITE);

            // file row
            DrawRectangle(0, ROW * 3, S_WIDTH*0.47, ROW,  Fade(LIGHTGRAY, 0.3f));
            DrawRectangle(S_WIDTH*0.47, ROW * 3, S_WIDTH, ROW,  Fade(BLUE, 0.3f));
            if (!file_is_dropped) {
                DrawText("Drop your files to this window!", 3, ROW*3 + ROW*0.5 - LARGE_TEXT*0.5, LARGE_TEXT, RED);
                DrawText("No file selected!", S_WIDTH*0.55, ROW*3 + ROW*0.5 - LARGE_TEXT*0.5, LARGE_TEXT, BLACK);
            }
            else if (multi_files_dropped) {
                DrawText("Only drop one file or folder!", 3, ROW*3 + ROW*0.5 - LARGE_TEXT*0.5, LARGE_TEXT, RED); 
                DrawText("No file selected!", S_WIDTH*0.55, ROW*3 + ROW*0.5 - LARGE_TEXT*0.5, LARGE_TEXT, BLACK);
            }
            else {
                DrawText("Dropped files:", 3, ROW*3 + ROW*0.5 - LARGE_TEXT*0.5, LARGE_TEXT, RED);
                DrawText(file_path, S_WIDTH*0.50, ROW*3 + ROW*0.5 - MEDIUM_TEXT*0.5, MEDIUM_TEXT, BLACK);
            }

            // run button 
            if (GuiButton((Rectangle){S_WIDTH - 200, S_HEIGHT - 200, 120, 30 }, "#191#Run Nanomux")) run_button = true;

            if (run_button)
            {
                int result = GuiMessageBox((Rectangle){ S_WIDTH / 2, S_HEIGHT / 2, 250, 100 },
                    "#191#Message Box", "Hi! This is a message!", "Nice;Cool");

                if (result >= 0) run_button = false;
            }

        EndDrawing();

    }
    CloseWindow();          
    return 0;
}