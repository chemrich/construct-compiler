from google import genai
import os

api_key = os.environ.get("GOOGLE_API_KEY")

if api_key:
    client = genai.Client(api_key=api_key)
    try:
        response = client.models.generate_content(
            model='gemini-2.5-flash',
            contents='Why is the sky blue?',
        )
        print(response.text)
    except Exception as e:
        print(f"API call failed: {e}")
else:
    print("GOOGLE_API_KEY not set.")
