%% Video Request %%

%  This function asks the user if he wants to create a video of the last
%  simulation. If the answer is 'yes' it creates it and stores it into the
%  videoDirectory defined in SetDirectories.

video_answer = input("\n\nCreate the video of the last simulation? \n 1 = yes \n 2 = no\n\nAnswer: ");

if video_answer == 1
    disp(["Video Processing ..."]);
    create_video;
elseif video_answer == 2
    disp(["Ok... Your simulation is complete"]);
else
    disp(["Error ...Close the application!"]);
end
